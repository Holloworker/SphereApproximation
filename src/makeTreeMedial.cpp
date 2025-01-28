/*************************************************************************\

                             C O P Y R I G H T

  Copyright 2003 Image Synthesis Group, Trinity College Dublin, Ireland.
                        All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following paragraphs appear in all copies.


                             D I S C L A I M E R

  IN NO EVENT SHALL TRININTY COLLEGE DUBLIN BE LIABLE TO ANY PARTY FOR
  DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING,
  BUT NOT LIMITED TO, LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE
  AND ITS DOCUMENTATION, EVEN IF TRINITY COLLEGE DUBLIN HAS BEEN ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGES.

  TRINITY COLLEGE DUBLIN DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE.  THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND TRINITY
  COLLEGE DUBLIN HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
  ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted at the following e-mail addresses:

          Gareth_Bradshaw@yahoo.co.uk    isg@cs.tcd.ie

  Further information about the ISG and it's project can be found at the ISG
  web site :

          isg.cs.tcd.ie

\**************************************************************************/

/*
    tests sphere trees generated with various methods
*/
#include <stdio.h>
#include <string.h>

#include "Surface/Surface.h"
#include "Surface/OBJLoader.h"
#include "API/MSGrid.h"
#include "MedialAxis/Voronoi3D.h"
#include "API/SEConvex.h"
#include "API/SESphPt.h"
#include "API/VFAdaptive.h"
#include "API/SSIsohedron.h"
#include "API/SRMerge.h"
#include "API/SRExpand.h"
#include "API/SRBurst.h"
#include "API/SRComposite.h"
#include "API/SFWhite.h"
#include "API/REMaxElim.h"
#include "API/STGGeneric.h"
#include "API/SOSimplex.h"
#include "API/SOBalance.h"
#include "EvalTree.h"
#include "VerifyModel.h"
#include "DecodeParam.h"
#include <iostream>
#include <fstream>
#include <Python.h>
#include "ManifoldPlus/src/Parser.h"
#include "ManifoldPlus/src/types.h"
#include <igl/writeOBJ.h>
#include "ManifoldPlus/src/IO.h"
#include "ManifoldPlus/src/Manifold.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/decimate.h>
#include <filesystem>
#include <mutex>
namespace fs = std::filesystem;
/*
    options and their default values
*/
int testerLevels =
    -1;         //  number of levels for NON-CONVEX, -1 uses CONVEX tester
int branch = 8; //  branching factor of the sphere-tree
int depth = 3;  //  depth of the sphere-tree
int numCoverPts = 5000; //  number of test points to put on surface for coverage
int minCoverPts = 5;    //  minimum number of points per triangle for coverage
int initSpheres = 500;  //  initial spheres in Voronoi
float erFact = 2;       //  error reduction factor for adaptive Voronoi
int spheresPerNode = 100; //  minimum number of spheres per node
float surface_num = 0.5;  //  percents of faces left(from 0 to 1)
bool verify = false;      //  verify model before construction
bool nopause = false;     //  will we pause before starting
bool eval = false;        //  do we evaluate the sphere-tree after construction
bool useMerge = false;    //  do we include the MERGE algorithm
bool useBurst = false;    //  do we include the BURST algorithm
bool useExpand = false;   //  do we include the EXPAND algorithm
bool Reduce_faces = false;

enum Optimisers { NONE, SIMPLEX, BALANCE };
int optimise = NONE;   //  which optimiser should we use
float balExcess = 0.0; //  % increase error allowed for BALANCE opt algorithm
int maxOptLevel =
    -1; //  maximum sphere-tree level to apply optimiser (0 does first set only)

#define DO_ALL_BELOW branch * 3
#define MAX_ITERS_FOR_VORONOI 1.50 * spheresPerNode;

/*
    info for decoding parameters
*/
const IntParam intParams[] = {{"testerLevels", &testerLevels},
                              {"maxOptLevel", &maxOptLevel},
                              {"branch", &branch},
                              {"depth", &depth},
                              {"numCover", &numCoverPts},
                              {"minCover", &minCoverPts},
                              {"initSpheres", &initSpheres},
                              {"minSpheres", &spheresPerNode},
                              {NULL, NULL}};

const FloatParam floatParams[] = {
    {"erFact", &erFact}, {"surface_num", &surface_num}, {"balExcess", &balExcess}, {NULL, NULL}};

const BoolParam boolParams[] = {{"verify", &verify, TRUE},
                                {"reduce_faces", &Reduce_faces, TRUE},
                                {"nopause", &nopause, TRUE},
                                {"eval", &eval, TRUE},
                                {"merge", &useMerge, TRUE},
                                {"burst", &useBurst, TRUE},
                                {"expand", &useExpand, TRUE},
                                {NULL, NULL}};

/*
    forward declarations
*/
void waitForKey();
int error(const char *errorMsg, const char *errorMsg1 = NULL);
void writeOptions(FILE *f);
void decodeOptions(int argc, const char *argv[]);
bool constructTree(const char *file);

void writeArrayToNewFile(const Array<Surface::Point>& arr, const std::string& filename) {
    // 创建并打开文件
    std::ofstream file(filename, std::ios::out | std::ios::trunc);
    
    if (!file.is_open()) {
        std::cerr << "Failed to create file: " << filename << std::endl;
        return;
    }

    // 将数组内容写入文件
    for (int i = 0; i < arr.getSize(); ++i) {
        file << arr.index(i).p.x;
        file << " ";
        file << arr.index(i).p.y;
        file << " ";
        file << arr.index(i).p.z;
        file << std::endl;
    }

    // 关闭文件
    file.close();
    if (file.fail()) {
        std::cerr << "Error writing to file: " << filename << std::endl;
    }
}
/*
    MAINLINE
*/
int main(int argc, const char *argv[]) {
  printf("MakeTreeMedial - Gareth Bradshaw Feb 2003\n");
  
  \
  /*
     parse command line
  */
  decodeOptions(argc, argv);
  if (!useMerge && !useBurst && !useExpand)
    useMerge = useExpand = TRUE; //  default to combined algorithm
  writeOptions(stdout);

  /*
     look for filenames and construct trees
  */
  int numFiles = 0;
  for (int i = 1; i < argc; i++) {
    if (argv[i] != NULL) {
      constructTree(argv[i]);
      numFiles++;
    }
  }

  /*
     check we had a file name
  */
  if (numFiles == 0)
    error("no files given :(");

  waitForKey();
}

/*
    construct sphere-tree for the model
*/
bool constructTree(const char *file) {
  std::string working_file(file);
  std::string manifold_file;
  std::string simplified_file;
  Py_Initialize();
  PyObject* pModule = PyImport_ImportModule("trimesh");
  if (!pModule) {
    PyErr_Print();
    std::cerr << "Failed to load Python module 'trimesh'" << std::endl;
    Py_Finalize();
    return false;
  }
  PyObject *args = PyTuple_New(1);
  PyObject *kwargs = PyDict_New();
  PyObject *load_func = PyObject_GetAttrString(pModule, "load"); 
  PyObject *file_path = PyUnicode_FromString(file);
  PyTuple_SetItem(args, 0, file_path);
  PyObject *process_key = PyUnicode_FromString("process");
  PyObject *process_value = PyBool_FromLong(0);  // Python's False
  PyDict_SetItem(kwargs, process_key, process_value);
  Py_DECREF(process_key);
  Py_DECREF(process_value);
  PyObject *mesh_or_scene = PyObject_Call(load_func, args, kwargs);
  PyObject* trimesh_class = PyObject_GetAttrString(pModule, "Trimesh");
  int is_instance = PyObject_IsInstance(mesh_or_scene, trimesh_class);
  int is_watertight = -1;
  if (is_instance == 1) {
    PyObject *is_watertight_obj = PyObject_GetAttrString(mesh_or_scene, "is_watertight");
    is_watertight = PyObject_IsTrue(is_watertight_obj);
  } else {
    PyObject *geometry = PyObject_GetAttrString(mesh_or_scene, "geometry");
    if (geometry && PyDict_Check(geometry)) {
        PyObject *key, *value;
        Py_ssize_t pos = 0;
        while (PyDict_Next(geometry, &pos, &key, &value)) {
          PyObject *mesh_name = PyObject_Str(key);
          const char *mesh_name_cstr = PyUnicode_AsUTF8(mesh_name);
          PyObject *is_watertight_obj = PyObject_GetAttrString(value, "is_watertight");
          is_watertight = PyObject_IsTrue(is_watertight_obj);
          Py_DECREF(is_watertight_obj);
          Py_DECREF(mesh_name);
          if (!is_watertight) {
            break;
          }
        }
    }
  }
  if (is_watertight == 1){
    printf("obj is watertight\n");
  }
  else if (is_watertight == 0){
    printf("obj is NOT watertight\n");
    printf("do manifold transition\n");
    MatrixD V, out_V;
	  MatrixI F, out_F;
	  ReadOBJ(working_file.c_str(), &V, &F);
    Manifold manifold;
    int dep = 8;
    manifold.ProcessManifold(V, F, dep, &out_V, &out_F);
    fs::path original_path(working_file);
    fs::path filename = original_path.filename();
    fs::path new_filename = "manifold_" + filename.string();
    fs::path new_path = original_path.parent_path() / new_filename;
    manifold_file = new_path.string();
    WriteOBJ(manifold_file.c_str(), out_V, out_F);
    while (!fs::exists(manifold_file)) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    working_file = manifold_file;
    std::cout << "num faces: " << out_F.rows() << std::endl;
  }
  else{
    printf("bad case!\n");
  }
  if (Reduce_faces){
    std::cout << "Simplifying the obj" << std::endl;
    std::cout << working_file << std::endl;
    MatrixD V_S;
	  MatrixI F_S;
    while (!fs::exists(working_file)) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    ReadOBJ(working_file.c_str(), &V_S, &F_S);
    Eigen::MatrixXd V_decimated;
    Eigen::MatrixXi F_decimated;
    Eigen::VectorXi J, I;
    int target_faces = F_S.rows() * surface_num;
    std::cout << "target faces: " << target_faces << std::endl;
    bool success = igl::decimate(V_S, F_S, target_faces, V_decimated, F_decimated, J, I);
    if (!success) {
      std::cerr << "Failed to reduce surfaces！" << std::endl;
      return -1;
    }
    fs::path original_path(working_file);
    fs::path filename = original_path.filename();
    fs::path new_filename = "simplified_" + filename.string();
    fs::path new_path = original_path.parent_path() / new_filename;
    simplified_file = new_path.string();
    igl::writeOBJ(simplified_file.c_str(), V_decimated, F_decimated);
    working_file = simplified_file;
  }
  Surface sur;
  std::cout << "filename in: " << working_file << std::endl;
        
  bool loaded = false;
  const char* file_str = working_file.c_str();
  if (strcompare("obj", file_str + strlen(file_str) - 3) == 0) {
    printf("Using OBJ loader\n");
    loaded = loadOBJ(&sur, file_str);
  } else {
    printf("Using SUR loader\n");
    loaded = sur.loadSurface(file_str);
  }
  if (!loaded) {
    throw std::runtime_error("Unable to load input file");
  }

  /*
      scale box
  */
  float cx_offset = 0.0;
  float cy_offset = 0.0;
  float cz_offset = 0.0;
  float boxScale = sur.fitIntoBox(1000, cx_offset, cy_offset, cz_offset);
  std::cout << "the offset is " << cx_offset << std::endl;
  std::cout << "the offset is " << cy_offset << std::endl;
  std::cout << "the offset is " << cz_offset << std::endl;
  /*
      make medial tester
  */
  MedialTester mt;
  mt.setSurface(sur);
  mt.useLargeCover = true;

  /*
      setup evaluator
  */
  SEConvex convEval;
  convEval.setTester(mt);
  SEBase *eval = &convEval;

  Array<Point3D> sphPts;
  SESphPt sphEval;
  if (testerLevels > 0) { //  <= 0 will use convex tester
    SSIsohedron::generateSamples(&sphPts, testerLevels - 1);
    sphEval.setup(mt, sphPts);
    eval = &sphEval;
    printf("Using concave tester (%d)\n\n", sphPts.getSize());
  }

  /*
      verify model
  */
  if (verify) {
    bool ok = verifyModel(sur);
    if (!ok) {
      printf("ERROR : model is not usable\n\n");
      return false;
    }
  }

  /*
      setup for the set of cover points
  */
  Array<Surface::Point> coverPts;
  MSGrid::generateSamples(&coverPts, numCoverPts, sur, TRUE, minCoverPts);
  writeArrayToNewFile(coverPts, "sample_points");
  printf("%d cover points\n", coverPts.getSize());

  /*
     Setup voronoi diagram
  */
  Point3D pC;
  pC.x = (sur.pMax.x + sur.pMin.x) / 2.0f;
  pC.y = (sur.pMax.y + sur.pMin.y) / 2.0f;
  pC.z = (sur.pMax.z + sur.pMin.z) / 2.0f;
  std::cout << "this maybe" <<std::endl;
  Voronoi3D vor;
  std::cout << "reach this maybe" <<std::endl;
  vor.initialise(pC, 1.5f * sur.pMin.distance(pC));

  /*
      setup adaptive Voronoi algorithm
  */
  VFAdaptive adaptive;
  adaptive.mt = &mt;
  adaptive.eval = eval;

  /*
      setup FITTER
  */
  SFWhite fitter;

  /*
      setup MERGE
  */
  SRMerge merger;
  merger.sphereEval = eval;
  merger.sphereFitter = &fitter;
  merger.useBeneficial = true;
  merger.doAllBelow = DO_ALL_BELOW;
  merger.setup(&vor, &mt);
  merger.vorAdapt = &adaptive;
  merger.initSpheres = initSpheres;
  merger.errorDecreaseFactor = erFact;
  merger.minSpheresPerNode = spheresPerNode;
  merger.maxItersForVoronoi = MAX_ITERS_FOR_VORONOI;

  /*
      setup BURST
  */
  SRBurst burster;
  burster.sphereEval = eval;
  burster.sphereFitter = &fitter;
  burster.useBeneficial = true;
  burster.doAllBelow = DO_ALL_BELOW;
  burster.setup(&vor, &mt);

  if (!useBurst) {
    //  setup adaptive algorithm
    burster.vorAdapt = &adaptive;
    burster.initSpheres = initSpheres;
    burster.errorDecreaseFactor = erFact;
    burster.minSpheresPerNode = spheresPerNode;
    burster.maxItersForVoronoi = MAX_ITERS_FOR_VORONOI;
  } else {
    // leave adaptive algorithm out as merge will do it for us
    burster.vorAdapt = NULL;
  }

  /*
      setup EXPAND generator
  */
  REMaxElim elimME;
  SRExpand expander;
  expander.redElim = &elimME;
  expander.setup(&vor, &mt);
  expander.errStep = 100;
  expander.useIterativeSelect = false;
  expander.relTol = 1E-5;
  if (!useMerge && !useBurst) {
    //  setup adaptive algorithm
    expander.vorAdapt = &adaptive;
    expander.initSpheres = initSpheres;
    expander.errorDecreaseFactor = erFact;
    expander.minSpheresPerNode = spheresPerNode;
    expander.maxItersForVoronoi = MAX_ITERS_FOR_VORONOI;
  } else {
    // leave adaptive algorithm out as previous algs will do it for us
    expander.vorAdapt = NULL;
  }

  /*
    setup the COMPOSITE algorithm
  */
  SRComposite composite;
  composite.eval = eval;
  if (useMerge)
    composite.addReducer(&merger);
  if (useBurst)
    composite.addReducer(&burster);
  if (useExpand)
    composite.addReducer(&expander);

  /*
      setup simplex optimiser in case we want it
  */
  SOSimplex simOpt;
  simOpt.sphereEval = eval;
  simOpt.sphereFitter = &fitter;

  /*
      setup balance optimiser to throw away spheres as long as
      increase in error is less than balExcess e.g. 0.05 is 1.05%
  */
  SOBalance balOpt;
  balOpt.sphereEval = eval;
  balOpt.optimiser = &simOpt;
  balOpt.V = 0.0f;
  balOpt.A = 1;
  balOpt.B = balExcess;

  /*
      setup SphereTree constructor - using dynamic construction
  */
  STGGeneric treegen;
  treegen.eval = eval;
  treegen.useRefit = true;
  treegen.setSamples(coverPts);
  treegen.reducer = &composite;
  treegen.maxOptLevel = maxOptLevel;
  if (optimise == SIMPLEX)
    treegen.optimiser = &simOpt;
  else if (optimise == BALANCE)
    treegen.optimiser = &balOpt;

  /*
      make sphere-tree
  */
  SphereTree tree;

  tree.setupTree(branch, depth + 1);

  waitForKey();
  treegen.constructTree(&tree);

  /*
     save sphere-tree
  */
  char sphereFile[1024];
  sprintf(sphereFile, "%s-medial.sph", working_file.c_str());
  if (tree.saveSphereTree(sphereFile, 1.0f / boxScale, branch, depth, cx_offset, cy_offset, cz_offset)) {
    Array<LevelEval> evals;
    if (eval) {
      evaluateTree(&evals, tree, eval);
      writeEvaluation(stdout, evals);
    }

    FILE *f = fopen(sphereFile, "a");
    if (f) {
      fprintf(f, "\n\n");
      writeOptions(f);
      fprintf(f, "\n\n");
      writeEvaluation(f, evals);
      fclose(f);
    }

    return true;
  } else {
    return false;
  }
}

/*
    error handler
*/
int error(const char *errorMsg, const char *errorMsg1) {
  if (errorMsg1)
    printf("ERROR : %s (%s)\n", errorMsg, errorMsg1);
  else
    printf("ERROR : %s\n", errorMsg);

  waitForKey();
  exit(-1);
}

void waitForKey() {
  if (!nopause) {
    printf("Press ENTER to continue....\n");
    char buffer[2];
    fread(buffer, 1, 1, stdin);
  }
}

/*
    command line parsing
*/
void decodeOptions(int argc, const char *argv[]) {
  //  decode the normal parameters
  decodeIntParam(argc, argv, intParams);
  decodeFloatParam(argc, argv, floatParams);
  decodeBoolParam(argc, argv, boolParams);

  //  decode extra parameters
  for (int i = 1; i < argc; i++) {
    if (argv[i] && strcompare(argv[i], "-optimise") == 0) {
      if (argv[i + 1] == NULL)
        error("decoding parameter", "-optimise");
      else if (strcompare(argv[i + 1], "off") == 0)
        optimise = NONE;
      else if (strcompare(argv[i + 1], "on") == 0 ||
               strcompare(argv[i + 1], "simplex") == 0)
        optimise = SIMPLEX;
      else if (strcompare(argv[i + 1], "balance") == 0)
        optimise = BALANCE;
      else
        error("decoding parameter", "-optimise");

      argv[i] = argv[i + 1] = NULL;
    }
  }
}

void writeOptions(FILE *f) {
  fprintf(f, "Options :\n");

  //  write out standard parameters
  writeParam(f, intParams);
  writeParam(f, floatParams);
  writeParam(f, boolParams);

  /*
    if (testerLevels <= 0)
      fprintf(f, "\ttesterLevels\t: Using CONVEX Tester\n");
    else{
      int numFaces = 20 * (int)pow(4, testerLevels-1);
      int numVerts = numFaces/2 + 2;
      fprintf(f, "\ttesterLevels\t: %d (%d points)\n", testerLevels, numVerts);
      }*/

  //  write out extra parameters
  if (optimise == SIMPLEX)
    fprintf(f, "\toptimise\t: SIMPLEX\n");
  else if (optimise == BALANCE)
    fprintf(f, "\toptimise\t: BALANCE (%.2f)\n", balExcess);
  else
    fprintf(f, "\toptimise\t: OFF\n");

  if (optimise != NONE)
    fprintf(f, "\tmaxOptLevel\t: %d\n", maxOptLevel);

  fprintf(f, "\n\n");
}
