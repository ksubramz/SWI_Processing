/**
 * @file swi_process.cpp
 * @brief Processing for Susceptibility Weighted Imaging Pipeline
 *
 */
/*
 * Original Author: Krish Subramaniam
 * CVS Revision Info:
 *    $Author$
 *    $Date$
 *
 * Copyright (C) 2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

//FreeSurfer C includes
extern "C"
{
#include "fsenv.h"
#include "mri.h"
#include "diag.h"
#include "DICOMRead.h"
};
#include <stdexcept>
#include <vector>
#include <iterator>
#include <algorithm>
#include "cmd_line_interface.h"
char *Progname;

// static function declarations
// forward declaration
struct IoParams;
using namespace std;

/*
  IO structure which store the command line args
*/
struct IoParams
{
  std::string magfile;
  std::string phasefile;
  std::string swioutput;
  float stddev;
  int phasemult;
  float phasecutoff;
  int miplevel;

  IoParams();
  void parse(int ac, char* av[]);
};

int main(int argc, char*argv[])
{
  IoParams params;
  try
  {
    params.parse(argc,argv);
  }
  catch (std::exception& excp)
  {
    std::cerr << "ERROR: Exception caught while parsing the command-line\n"
    << excp.what() << std::endl;
    exit(1);
  }

  std::cout << "Reading the input file(s)\n";
  MRI *mrimag=NULL, *mriphase=NULL;
  
  mrimag = MRIread(params.magfile.c_str() );
  std::cout << "Read the files..\n";
  if ( mrimag == NULL )
  {
    std::cerr << "ERROR: The magnitude image can't be read. Check the file\n";
    exit(1);
  }
  mriphase = MRIread(params.phasefile.c_str() );
  if ( mriphase == NULL )
  {
    std::cerr << "ERROR: The phase image can't be read. Check the file\n";
    exit(1);
  }

  float mval, pval, outval, _ov =0.0, v;
  float vmin, vmax;
  MRI *outimg = MRIclone(mriphase, NULL);
  // the following are the various intermediary images
  MRI *_hpimg = MRIclone(mriphase, NULL);
  MRI *_tmpimg = MRIclone(mriphase, NULL);

  std::cout << "Performing Gaussian smoothing with stddev " << params.stddev << " ..\n";
  MRIgaussianSmooth(mriphase, params.stddev, 1, _tmpimg);
  //MRIwrite(_tmpimg, "gaussian.mgz");

  std::cout << "Subtracting the smoothed image from the phase image..\n";
  MRIsubtract(mriphase, _tmpimg, _tmpimg);
  /*MRI_fft_highpass(mriphase, _tmpimg, 20);*/
  //MRIwrite(_tmpimg, "highpass.mgz");
  //
  MRIvalRange(_tmpimg, &vmin, &vmax);
  std::cout << "Min and max of the highpass image " << vmin << " " << vmax << "\n";

  if ( params.phasecutoff < vmin )
    params.phasecutoff = vmin;
  std::cout << "Performing phase mask multiplication with cutoff " << params.phasecutoff << " ..\n";
  std::cout << "Number of phase multiplications " << params.phasemult << " ..\n";
  for ( int i=0; i < _hpimg->width; i++)
    for ( int j=0; j < _hpimg->height; j++)
      for ( int k=0; k < _hpimg->depth; k++)
      {
        mval = MRIgetVoxVal(mrimag, i, j, k, 0);
        pval = MRIgetVoxVal(_tmpimg, i, j, k, 0);
        // do phase cut off
        if ( pval > 0.0 )
          _ov = 1.0;
        else if ( pval <= params.phasecutoff )
          _ov = 0.0;
        else
          _ov = 1.0 - (pval / params.phasecutoff );
        
        // do phase multiplications
        outval=1.0;
        for (int phm=0; phm < params.phasemult; phm++) 
          outval = outval * _ov;
        v = outval * mval;
        MRIsetVoxVal(_hpimg, i, j, k, 0, v);
      }

  // minimum Intensity Projection along y
  std::cout << "Performing Minimum Intensity Projection along y direction with levels=" << params.miplevel << " ..\n";
  int tmp_idx;
  std::vector<float> vals(params.miplevel);
  float mip_val;
  for ( int i=0; i < _hpimg->width; i++)
    for ( int j=0; j < _hpimg->height; j++)
      for ( int k=0; k < _hpimg->depth; k++)
        {
          for ( int jc=j; jc < j+params.miplevel; jc++)
          {
            tmp_idx = jc;
            // bounds check
            if ( jc < 0 )
              tmp_idx = 0;
            if ( jc >= _hpimg->height-1 )
              tmp_idx = _hpimg->height - 1;
            vals.push_back( MRIgetVoxVal(_hpimg, i, tmp_idx, k, 0));
          }

          mip_val = *(std::min_element( vals.begin(), vals.end() ));
          vals.clear();
          MRIsetVoxVal(outimg, i, j, k, 0, mip_val);
	}
  std::cout << "Writing the swi processed image\n";
  MRIwrite(outimg, params.swioutput.c_str() );
  // Freeing all
  MRIfree(&_hpimg);
  MRIfree(&mrimag);
  MRIfree(&mriphase);
  MRIfree(&outimg);
  MRIfree(&_tmpimg);
}

// These are the default values for the command line arguments 
IoParams::IoParams()
{
  stddev      = 2.0;
  phasecutoff = -std::numeric_limits<float>::max();
  phasemult   = 4;
  miplevel    = 4;
}

void IoParams::parse(int argc, char* argv[])
{
  CCmdLineInterface cmd(argv[0]);

  cmd.AddOptionString("mag_file", &magfile, "The magnitude image ( Output from the PRELUDE program)");
  cmd.AddOptionString("phase_file", &phasefile, "The phase image ( Output from the PRELUDE program)");
  cmd.AddOptionString("swi_output", &swioutput, "Name of the SWI processed output image");
  cmd.AddOptionFloat("stddev", &stddev, "Specify the standard deviation of the Gaussian Smoothing Filter. Default is 2.0");
  cmd.AddOptionInt("phase_multiplications", &phasemult, "Specify the number of phase multiplications. Default is 4");
  cmd.AddOptionFloat("phase_mask_cutoff", &phasecutoff, "Specify the negative phase mask cutoff frequency ( in radianѕ). Default is the minimum value of the phase image.");
  cmd.AddOptionInt("mip_level", &miplevel, "Specify the number of levels of mIP across the y direction. Default is 4");
  if ( argc == 1 )
  {
    std::cout << "\n"
    " Process the Susceptibility-weighted images. Make sure the inputs to this program is after the phase unwrapping step using PRELUDE\n";
    cmd.PrintHelp();
    exit(0);
  }

  cmd.Parse(argc, argv);
  
}


