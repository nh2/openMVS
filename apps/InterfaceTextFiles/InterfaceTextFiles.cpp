/*
 * InterfaceTextFiles.cpp
 *
 * Copyright (c) 2014-2015 SEACAVE
 *
 * Author(s):
 *
 *      cDc <cdc.seacave@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Additional Terms:
 *
 *      You are required to preserve legal notices and author attributions in
 *      that material or in the Appropriate Legal Notices displayed by works
 *      containing it.
 */

#include "../../libs/MVS/Common.h"
#include "../../libs/MVS/Scene.h"
#define LOG_OUT() GET_LOG()
#define LOG_ERR() GET_LOG()
#include "Util.h"
#include <map>
#include <boost/program_options.hpp>


// D E F I N E S ///////////////////////////////////////////////////

#define APPNAME _T("InterfaceTextFiles")
#define MVS_EXT _T(".mvs")
#define VSFM_EXT _T(".nvm")


// S T R U C T S ///////////////////////////////////////////////////

struct Point3D {
  Point3D() : point3DId(-1), x(0), y(0), z(0), r(0), g(0), b(0), error(0) {}

  int point3DId;
  double x, y, z;
  int r, g, b;
  double error;

  std::vector<uint32_t> imageIndices;
};

namespace OPT {
String strInputFileName;
String strUndistortedImagesDir;
String strOutputFileName;
String strOutputImageFolder;
unsigned nArchiveType;
int nProcessPriority;
unsigned nMaxThreads;
String strConfigFileName;
boost::program_options::variables_map vm;
} // namespace OPT

// initialize and parse the command line parameters
bool Initialize(size_t argc, LPCTSTR* argv)
{
  // initialize log and console
  OPEN_LOG();
  OPEN_LOGCONSOLE();

  // group of options allowed only on command line
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,h", "produce this help message")
    ("working-folder,w", boost::program_options::value<std::string>(&WORKING_FOLDER), "working directory (default current directory)")
    ("config-file,c", boost::program_options::value<std::string>(&OPT::strConfigFileName)->default_value(APPNAME _T(".cfg")), "file name containing program options")
    ("archive-type", boost::program_options::value<unsigned>(&OPT::nArchiveType)->default_value(2), "project archive type: 0-text, 1-binary, 2-compressed binary")
    ("process-priority", boost::program_options::value<int>(&OPT::nProcessPriority)->default_value(-1), "process priority (below normal by default)")
    ("max-threads", boost::program_options::value<unsigned>(&OPT::nMaxThreads)->default_value(0), "maximum number of threads (0 for using all available cores)")
    #if TD_VERBOSE != TD_VERBOSE_OFF
    ("verbosity,v", boost::program_options::value<int>(&g_nVerbosityLevel)->default_value(
      #if TD_VERBOSE == TD_VERBOSE_DEBUG
      3
      #else
      2
      #endif
      ), "verbosity level")
    #endif
    ;

  // group of options allowed both on command line and in config file
  boost::program_options::options_description config("Main options");
  config.add_options()
    ("input-file,i", boost::program_options::value<std::string>(&OPT::strInputFileName), "input filename containing camera poses and image list")
    ("undistorted-images-dir,u", boost::program_options::value<std::string>(&OPT::strUndistortedImagesDir), "directory containing undistorted images, as exported by colmap")
    ("output-file,o", boost::program_options::value<std::string>(&OPT::strOutputFileName), "output filename for storing the mesh")
    ("output-image-folder", boost::program_options::value<std::string>(&OPT::strOutputImageFolder)->default_value("undistorted_images"), "output folder to store undistorted images")
    ;

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(config);

  boost::program_options::options_description config_file_options;
  config_file_options.add(config);

  boost::program_options::positional_options_description p;
  p.add("input-file", -1);

  try {
    // parse command line options
    boost::program_options::store(boost::program_options::command_line_parser((int)argc, argv).options(cmdline_options).positional(p).run(), OPT::vm);
    boost::program_options::notify(OPT::vm);
    INIT_WORKING_FOLDER;
    // parse configuration file
    std::ifstream ifs(MAKE_PATH_SAFE(OPT::strConfigFileName));
    if (ifs) {
      boost::program_options::store(parse_config_file(ifs, config_file_options), OPT::vm);
      boost::program_options::notify(OPT::vm);
    }
  }
  catch (const std::exception& e) {
    LOG(e.what());
    return false;
  }

  // initialize the log file
  OPEN_LOGFILE(MAKE_PATH(APPNAME _T("-")+Util::getUniqueName(0)+_T(".log")).c_str());

  // print application details: version and command line
  Util::LogBuild();
  LOG(_T("Command line:%s"), Util::CommandLineToString(argc, argv).c_str());

  // validate input
  Util::ensureValidPath(OPT::strInputFileName);
  Util::ensureUnifySlash(OPT::strInputFileName);
  Util::ensureValidPath(OPT::strUndistortedImagesDir);
  Util::ensureUnifySlash(OPT::strUndistortedImagesDir);
  if (OPT::vm.count("help") || OPT::strInputFileName.IsEmpty() || OPT::strUndistortedImagesDir.IsEmpty()) {
    boost::program_options::options_description visible("Available options");
    visible.add(generic).add(config);
    GET_LOG() << visible;
  }
  if (OPT::strInputFileName.IsEmpty())
    return false;
  if (OPT::strUndistortedImagesDir.IsEmpty())
    return false;

  // initialize optional options
  Util::ensureValidPath(OPT::strOutputFileName);
  Util::ensureUnifySlash(OPT::strOutputFileName);
  Util::ensureUnifySlash(OPT::strOutputImageFolder);
  Util::ensureDirectorySlash(OPT::strOutputImageFolder);
  if (OPT::strOutputFileName.IsEmpty())
    OPT::strOutputFileName = Util::getFullFileName(OPT::strInputFileName) + MVS_EXT;

  // initialize global options
  Process::setCurrentProcessPriority((Process::Priority)OPT::nProcessPriority);
  #ifdef _USE_OPENMP
  if (OPT::nMaxThreads != 0)
    omp_set_num_threads(OPT::nMaxThreads);
  #endif

  #ifdef _USE_BREAKPAD
  // start memory dumper
  MiniDumper::Create(APPNAME, WORKING_FOLDER);
  #endif
  return true;
}

// finalize application instance
void Finalize()
{
  #if TD_VERBOSE != TD_VERBOSE_OFF
  // print memory statistics
  Util::LogMemoryInfo();
  #endif

  CLOSE_LOGFILE();
  CLOSE_LOGCONSOLE();
  CLOSE_LOG();
}

#define PBA_PRECISION float

namespace PBA {
template<class FT> struct CameraT_;
typedef CameraT_<PBA_PRECISION> Camera;
template<class FT> struct Point3D_;
typedef Point3D_<PBA_PRECISION> Point3D;
} // namespace PBA

namespace MVS {
} // namespace MVS

int main(int argc, LPCTSTR* argv)
{
  std::cerr << "main" << std::endl;

  #ifdef _DEBUGINFO
  // set _crtBreakAlloc index to stop in <dbgheap.c> at allocation
  _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);// | _CRTDBG_CHECK_ALWAYS_DF);
  #endif

  if (!Initialize(argc, argv))
    return EXIT_FAILURE;

  TD_TIMER_START();

  // read VisualSFM input data
  // std::vector<PBA::Camera> cameras;
  // std::vector<PBA::Point3D> vertices;
  // std::vector<PBA::Point2D> measurements; // the array of 2D projections (only inliers)
  // std::vector<int> correspondingPoint; // 3D point index corresponding to each 2D projection
  // std::vector<int> correspondingView; // and camera index
  // std::vector<std::string> names;
  // std::vector<int> pointColors;

  // Have to fill:
  //   - cameras
  //   - vertices
  //   - measurements
  //   - correspondingPoint
  //   - correspondingView
  //   - names
  //   - pointColors

  std::vector<std::string> imageFileNames;

  std::ifstream in(MAKE_PATH_SAFE(OPT::strInputFileName + "/images.txt"));

  // Parse file line by line
  std::string line;
  while (std::getline(in, line))
  {
    // Skip comment lines (the headers at the top are also comments)
    if (line.length() >= 1 && line[0] == '#') continue;

    // Colmap images.txt files have the following format:
    //    # Image list with two lines of data per image:
    //    #   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME
    //    #   POINTS2D[] as (X, Y, POINT3D_ID)

    int imageId;
    double qw, qx, qy, qz, tx, ty, tz;
    int cameraId;
    std::string name;

    // Read a
    //   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME
    // line; example:
    //   16 0.791455 0.0325857 -0.541978 -0.280709 2.18149 -2.30782 3.82575 1 IMG_7147.JPG
    {
      std::istringstream lineStream(line);
      if (!(lineStream >> imageId)) { throw std::runtime_error("expected IMAGE_ID"); }
      if (!(lineStream >> qw >> qx >> qy >> qz)) { throw std::runtime_error("expected image rotation quaternions QW, QX, QY, QZ"); }
      if (!(lineStream >> tx >> ty >> tz)) { throw std::runtime_error("expected image position vector TX, TY, TZ"); }
      if (!(lineStream >> cameraId)) { throw std::runtime_error("expected CAMERA_ID"); }
      if (!(lineStream >> name)) { throw std::runtime_error("expected image NAME"); }
      if (!lineStream.eof()) { throw std::runtime_error("expected end of line"); }

      // std::cout << "Image " << name << " ID = " << imageId << " Camera ID = " << cameraId
      //           << " rotation quaternion ( " << qw << " " << qx << " " << qy << " " << qz << " )"
      //           << " position ( " << tx << " " << ty << " " << tz << " )"
      //           << std::endl;
    }

    // Skip comment lines
    while (std::getline(in, line) && line.length() >= 1 && line[0] == '#') {};

    // Read a
    //   POINTS2D[] as (X, Y, POINT3D_ID)
    // line; example:
    //   1762.85 47.6561 -1 1936.35 72.2963 -1 2398.05 2412.79 8
    // where `-1` as POINT3D_ID means that this point is only a 2D feature point
    // and its 3D location could not be determined.
    {
      std::istringstream lineStream(line);
      while(!lineStream.eof()) {
        double x, y;
        int point3DId = -1;
        if (!(lineStream >> x >> y)) { throw std::runtime_error("expected X, Y"); }
        if (!(lineStream >> point3DId)) { throw std::runtime_error("expected POINT3D_ID"); }

        // std::cout << "Point3D ID = " << point3DId << " at X = " << x << " Y = " << y << std::endl;
      }
    }

    imageFileNames.push_back(name);
  }

  std::map<std::string, Matrix3x4> projMatrices;

  for (std::string name : imageFileNames) {

    // Read 3x4 projection matrix from undistortedImagesDir/NAME.proj_matrix.txt
    std::ifstream in(MAKE_PATH_SAFE(OPT::strUndistortedImagesDir + "/" + name + ".proj_matrix.txt"));

    Matrix3x4 projMatrix;
    in >> projMatrix(0, 0);
    in >> projMatrix(0, 1);
    in >> projMatrix(0, 2);
    in >> projMatrix(0, 3);
    in >> projMatrix(1, 0);
    in >> projMatrix(1, 1);
    in >> projMatrix(1, 2);
    in >> projMatrix(1, 3);
    in >> projMatrix(2, 0);
    in >> projMatrix(2, 1);
    in >> projMatrix(2, 2);
    in >> projMatrix(2, 3);

    // std::cout << "projMatrix for " << name << ":" << std::endl
    //           << projMatrix << std::endl;

    projMatrices.insert(std::pair<std::string, Matrix3x4>(name, projMatrix));
  }


  size_t numCameras = imageFileNames.size();

  MVS::Scene scene(OPT::nMaxThreads);
  scene.platforms.Reserve(numCameras);
  scene.images.Reserve(numCameras);
  scene.nCalibratedImages = 0;
  for (size_t idx=0; idx<numCameras; ++idx) {
    MVS::Image& image = scene.images.AddEmpty();
    std::string filename = imageFileNames[idx];
    image.name = MAKE_PATH_SAFE(OPT::strUndistortedImagesDir + "/" + filename);
    Util::ensureUnifySlash(image.name);
    image.name = MAKE_PATH_FULL(WORKING_FOLDER_FULL, image.name);
    std::cerr << image.name << std::endl;
    if (!image.ReloadImage(0, false)) {
      LOG("error: can not read image %s", image.name.c_str());
      return EXIT_FAILURE;
    }
    // set camera
    image.platformID = scene.platforms.GetSize();
    MVS::Platform& platform = scene.platforms.AddEmpty();
    MVS::Platform::Camera& camera = platform.cameras.AddEmpty();
    image.cameraID = 0;

    std::map<std::string, Matrix3x4>::iterator projMatrixIter = projMatrices.find(filename);
    assert(projMatrixIter != projMatrices.end());
    Matrix3x4 projMatrix = (*projMatrixIter).second;

    KMatrix K;
    RMatrix R;
    CMatrix C;
    MVS::DecomposeProjectionMatrix(projMatrix, K, R, C);

    double focalLengthX = K(0, 0);
    double focalLengthY = K(1, 1);
    camera.K = MVS::Platform::Camera::ComposeK<REAL,REAL>(focalLengthX, focalLengthY, image.width, image.height);
    camera.R = RMatrix::IDENTITY;
    camera.C = CMatrix::ZERO;
    // normalize camera intrinsics
    const REAL fScale(REAL(1)/MVS::Camera::GetNormalizationScale(image.width, image.height));
    camera.K(0, 0) *= fScale;
    camera.K(1, 1) *= fScale;
    camera.K(0, 2) *= fScale;
    camera.K(1, 2) *= fScale;
    // set pose
    image.poseID = platform.poses.GetSize();
    MVS::Platform::Pose& pose = platform.poses.AddEmpty();
    pose.R = R;
    pose.C = C;
    image.UpdateCamera(scene.platforms);
    ++scene.nCalibratedImages;
  }

  std::ifstream pointsIn(MAKE_PATH_SAFE(OPT::strInputFileName + "/points3D.txt"));

  std::vector<Point3D> points3d;

  while(std::getline(pointsIn, line)) {

    // Skip comment lines (the headers at the top are also comments)
    if (line.length() >= 1 && line[0] == '#') continue;

    // Colmap points3D.txt files have the following format:
    //    # 3D point list with one line of data per point:
    //    #   POINT3D_ID, X, Y, Z, R, G, B, ERROR, TRACK[] as (IMAGE_ID, POINT2D_IDX)

    Point3D p;

    // Read a
    //   POINT3D_ID, X, Y, Z, R, G, B, ERROR, TRACK[]
    // line; example:
    //   2661 1.59224 0.943826 3.27984 134 42 33 1.86371 2 1946 1 1993 3 2394 12 1794 17 631

    std::istringstream lineStream(line);
    if (!(lineStream >> p.point3DId)) { throw std::runtime_error("expected POINT3D_ID"); }
    if (!(lineStream >> p.x >> p.y >> p.z)) { throw std::runtime_error("expected X, Y, Z"); }
    if (!(lineStream >> p.r >> p.g >> p.b)) { throw std::runtime_error("expected R, G, B"); }
    // TODO check that r, g, b are within uint8
    if (!(lineStream >> p.error)) { throw std::runtime_error("expected ERROR"); }

    while (!lineStream.eof()) {
      uint32_t imageId;
      uint32_t point2DIndex;
      if (!(lineStream >> imageId)) { throw std::runtime_error("expected IMAGE_ID"); }
      if (!(lineStream >> point2DIndex)) { throw std::runtime_error("expected POINT2D_IDX"); }
      p.imageIndices.push_back(imageId);
    }

    points3d.push_back(p);
  }

  size_t numPoints3d = points3d.size();
  // std::cerr << "numPoints3d " << numPoints3d << std::endl;;
  scene.pointcloud.points.Reserve(numPoints3d);
  scene.pointcloud.colors.Reserve(numPoints3d * 3);
  scene.pointcloud.pointViews.Resize(numPoints3d);

  std::cerr << "images " << scene.images.size() << " " << numCameras << std::endl;
  for (size_t idx=0; idx<numPoints3d; ++idx) {
    Point3D p = points3d[idx];
    // std::cerr << "adding point " << p.x << " " << p.y << " " << p.z << std::endl;
    scene.pointcloud.points.AddConstruct(p.x, p.y, p.z);
    scene.pointcloud.colors.AddConstruct(p.r, p.g, p.b); // TODO check RGB order (or BGR?)

    MVS::PointCloud::ViewArr& views = scene.pointcloud.pointViews[idx];
    assert(p.imageIndices.size() > 0);
    for (uint32_t i : p.imageIndices) {
      // std::cerr << i << " ";
      views.InsertSort(i - 1);
    }
  }
  std::cerr << std::endl;

  // write OpenMVS input data
  scene.Save(MAKE_PATH_SAFE(OPT::strOutputFileName), (ARCHIVE_TYPE)OPT::nArchiveType);

  std::cerr << "Saving done" << std::endl;

  // VERBOSE("Input data imported: %u cameras, %u poses, %u images, %u vertices (%s)", numCameras, cameras.size(), cameras.size(), vertices.size(), TD_TIMER_GET_FMT().c_str());
  std::cerr << "Input data imported: " << scene.platforms.size() << " platforms" << std::endl;

  // if (!scene.Load(MAKE_PATH_SAFE(OPT::strOutputFileName))) {
  //   std::cerr << "loading scene failed" << std::endl;
  //   return 1;
  // }
  // std::cerr << "images:" << std::endl;
  // for (auto img : scene.images) {
  //   std::cerr << "- name: " << img.name << std::endl;
  // }
  // if (scene.pointcloud.IsEmpty()) {
  //   std::cerr << "error: empty initial point-cloud" << std::endl;
  // }
  // std::cerr << "pointcloud.points:" << std::endl;
  // for (auto p : scene.pointcloud.points) {
  //   std::cerr << "- " << p.x << " " << p.y << " " << p.z << std::endl;
  // }

  Finalize();
  return EXIT_SUCCESS;
}
/*----------------------------------------------------------------*/
