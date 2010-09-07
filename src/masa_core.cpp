// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Author$
// $Id$
//
// masa_core.cpp: this is the core set of routines -- the only functions that 
//                should be called by users -- all other files are internal
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <masa_internal.h>

using namespace MASA;

map<string, manufactured_solution*> masa_master_list; // global map between unique name and manufactured class
manufactured_solution* masa_master_pointer;           // pointer to currently selected manufactured solution

//
//  this function selects an already initialized manufactured class
//
int MASA::masa_select_mms(string name)
{
  string nametemp;
  int selector;
  // lets run though the list to check the variable does exist
  map<string,manufactured_solution*>::iterator it;
  it=masa_master_list.find(name);
  if(it != masa_master_list.end()) // found a name
    { 
      cout << "selected " << name << endl;
      masa_master_pointer=masa_master_list[name]; // set pointer to currently selected solution      
    }      
  else 
    {
      cout << "\nMASA ERROR: No such manufactured solution (" << name << ") has been initialized.\n";
      exit(1);
    } 

  return 0;
  
}// done with masa_select

// helper function (deprecated!)
int masa_v2o(void* obid, manufactured_solution** manfac)
{ 
  string name;
  manufactured_solution* temp;
  temp=static_cast<manufactured_solution*>(obid); // cast the void ptr back as an obj
  *manfac=temp;
  return 0;

}

// get list mms
//
// this function returns a vector of pointers to 
// manufactured solutions available
//
// adding function because it is called in several places

int get_list_mms(vector<manufactured_solution*>* anim)
{
  anim->push_back(new MASA_Test()); // test function
  
  // register solutions here
  anim->push_back(new heateq_1d_steady_const());
  anim->push_back(new heateq_2d_steady_const());
  anim->push_back(new heateq_3d_steady_const());

  anim->push_back(new heateq_1d_unsteady_const());
  anim->push_back(new heateq_2d_unsteady_const());
  anim->push_back(new heateq_3d_unsteady_const());

  anim->push_back(new heateq_1d_unsteady_var());
  anim->push_back(new heateq_2d_unsteady_var());
  anim->push_back(new heateq_3d_unsteady_var());

  anim->push_back(new heateq_1d_steady_var());
  anim->push_back(new heateq_2d_steady_var());
  anim->push_back(new heateq_3d_steady_var());

  anim->push_back(new euler_1d());
  anim->push_back(new euler_2d());
  anim->push_back(new euler_3d());

  anim->push_back(new sod_1d());

  anim->push_back(new navierstokes_2d_compressible());
  anim->push_back(new navierstokes_3d_compressible());

  anim->push_back(new axi_euler());
  anim->push_back(new axi_cns());

  return 0;

}

//
//  this function will initiate a masa manufactured class
//
int MASA::masa_init(string unique_name, string str)
{
  int flag=0;
  string name,temp;
  int error=1;
  vector<manufactured_solution*> anim;

  get_list_mms(&anim); //construct list 
  
  temp=str;
  masa_map(&temp); // remove lowercase, etc.

  for (vector<manufactured_solution*>::const_iterator it = anim.begin(); it != anim.end(); ++it) 
    {
      (*it)->return_name(&name); // get name
      error=temp.rfind(name);   // look for name -- must be identical to full name, after masa_map edits
      if (error!=string::npos) // found a value
	{
	  masa_master_list[unique_name]=*it; // write down name 
	  masa_master_pointer=*it; // set as currently active manufactured solution
	  flag=1;                 // set exit flag
	}
      else // strings not identical
	{
	  delete *it; // this calls the deconstructor
	}
    }// done with for loop 
  
  
  if(flag==1)
    {      
      //cout << "\nMASA got it\n";
    }
  else
    {
      cout << "\nMASA FATAL ERROR: No Manufactured Solutions of that Type\n";
      exit(1); // error code, terminate
    }
  
  return 0; // steady as she goes

}

//
//  this function will initiate a masa manufactured class
//
int MASA::masa_curr_mms(string* str)
{

  // lets run though the list to check the variable does exist
  masa_master_pointer->return_name(str);
  // cout << masa_master_list[masa_master_pointer] << endl;    

  return 0;
}

//
//  this function will initiate a masa manufactured class
//
int MASA::masa_list_mms()
{
  string str;

  // output the size of the map
  cout << "Number of initialized solutions: " << masa_master_list.size() << endl;
  for(map<string,manufactured_solution*>::iterator iter = masa_master_list.begin(); iter != masa_master_list.end(); iter++)
    {
      (iter->second)->return_name(&str);
      cout << iter->first << " : " << str << endl;
    }
  return 0;
}

//
// function that searches all registered masa solutions
// for a selected manufactured solution
// (deprecated)
int MASA::masa_getid(void** objid,string str)
{
  int flag=0;
  string name;
  vector<manufactured_solution*> anim;

  get_list_mms(&anim); //construct list 
 
  // masa_map();  

  for (vector<manufactured_solution*>::const_iterator it = anim.begin(); it != anim.end(); ++it) {
    (*it)->return_name(&name); // get name
    if(str.compare(name) == 0) // strings are identical
      {
	*objid=(*it); // cast object pointer as void
	flag=1;      // set flag
      }
    else // strings not identical
      delete *it; // this calls the deconstructor
    
  }// done with for loop 

  if(flag==1)
    {      
      //cout << "\nMASA got it\n";
    }
  else
    {
      cout << "\nMASA FATAL ERROR: No Manufactured Solutions of that Type\n";
      exit(1); // error code, terminate
    }

  return 0; // steady as she goes
}// done with masa getid

//
// function that prints all registered masa solutions
//
int MASA::masa_printid()
{
  vector<manufactured_solution*> anim;
  string name;

  get_list_mms(&anim); //construct list 

  cout << endl;
  // masa_map();
  for (vector<manufactured_solution*>::const_iterator it = anim.begin(); it != anim.end(); ++it) {
    (*it)->return_name(&name); // get name
    cout << endl << name;
    delete *it; // this calls the deconstructor
  }// done with for loop 
  
  return 0; // steady as she goes
}// done with masa print id

int MASA::masa_set_param(string param,double paramval)
{
  masa_master_pointer->set_var(param,paramval);
  return 0;
}

//
// Set all parameters to default values
//
int MASA::masa_init_param()
{
  masa_master_pointer->init_var();
  return 0;
}

//
// Function that returns value of parameter selected by string
// 

int MASA::masa_get_param(string param,double *paramval)
{

  masa_master_pointer->get_var(param,paramval);

  return 0;
}

int MASA::masa_display_param()
{

  masa_master_pointer->display_var();

  return 0;
}

/* ------------------------------------------------
 *
 *         1D functions
 *
 * -----------------------------------------------
 */ 


  // --------------------------------
  // source terms
  // --------------------------------

int MASA::masa_eval_t_source(double x,double* field) //x 
{
  *field=masa_master_pointer->eval_q_t(x);
  return 0;
}

int MASA::masa_eval_t_source(double x,double t,double* field) //x,t
{
  *field=masa_master_pointer->eval_q_t(x,t);
  return 0;
}

int MASA::masa_eval_u_source(double x,double* field)
{
  *field=masa_master_pointer->eval_q_u(x);
  return 0;
}

int MASA::masa_eval_rho_source(double x,double* field)
{
  *field=masa_master_pointer->eval_q_rho(x);
  return 0;
}

int MASA::masa_eval_rho_u_source(double x,double t,double* field)
{
  *field=masa_master_pointer->eval_q_rho_u(x,t);
  return 0;
}

int MASA::masa_eval_e_source(double x,double* field)
{
  *field=masa_master_pointer->eval_q_e(x);
  return 0;
}

  // --------------------------------
  // analytical terms
  // --------------------------------

int MASA::masa_eval_t_an(double x,double* field)
{
  *field=masa_master_pointer->eval_an_t(x);
  return 0;
}

int MASA::masa_eval_u_an(double x,double* field)
{
  *field=masa_master_pointer->eval_an_u(x);
  return 0;
}

int MASA::masa_eval_p_an(double x,double* field)
{
  *field=masa_master_pointer->eval_an_p(x);
  return 0;
}

int MASA::masa_eval_rho_an(double x,double* field)
{
  *field=masa_master_pointer->eval_an_rho(x);
  return 0;
}


/* ------------------------------------------------
 *
 *         2D functions
 *
 * -----------------------------------------------
 */ 

  // --------------------------------
  // source terms
  // --------------------------------

int MASA::masa_eval_t_source(double x,double y,double t,double* field)
{
  *field=masa_master_pointer->eval_q_t(x,y,t);
  return 0;
}

int MASA::masa_eval_u_source(double x,double y,double* field)
{

  *field=masa_master_pointer->eval_q_u(x,y);
  return 0;
}

int MASA::masa_eval_v_source(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_q_v(x,y);
  return 0;
}

int MASA::masa_eval_w_source(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_q_w(x,y);
  return 0;
}

int MASA::masa_eval_rho_source(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_q_rho(x,y);
  return 0;
}

int MASA::masa_eval_e_source(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_q_e(x,y);
  return 0;
}

  // --------------------------------
  // analytical terms
  // --------------------------------

int MASA::masa_eval_t_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_t(x,y);
  return 0;
}

int MASA::masa_eval_u_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_u(x,y);
  return 0;
}

int MASA::masa_eval_v_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_v(x,y);
  return 0;
}

int MASA::masa_eval_w_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_w(x,y);
  return 0;
}

int MASA::masa_eval_p_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_p(x,y);
  return 0;
}

int MASA::masa_eval_rho_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_rho(x,y);
  return 0;
}


/* ------------------------------------------------
 *
 *         3D functions
 *
 * -----------------------------------------------
 */ 

  // --------------------------------
  // source terms
  // --------------------------------

int MASA::masa_eval_t_source(double x,double y,double z,double t,double* field)
{
  *field=masa_master_pointer->eval_q_t(x,y,z,t);
  return 0;
}

int MASA::masa_eval_u_source(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_q_u(x,y,z);
  return 0;
}

int MASA::masa_eval_v_source(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_q_v(x,y,z);
  return 0;
}

int MASA::masa_eval_w_source(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_q_w(x,y,z);
  return 0;
}

int MASA::masa_eval_rho_source(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_q_rho(x,y,z);
  return 0;
}

int MASA::masa_eval_e_source(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_q_e(x,y,z);
  return 0;
}

  // --------------------------------
  // analytical terms
  // --------------------------------

int MASA::masa_eval_t_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_t(x,y,z);
  return 0;
}

int MASA::masa_eval_u_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_u(x,y,z);
  return 0;
}

int MASA::masa_eval_v_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_v(x,y,z);
  return 0;
}

int MASA::masa_eval_w_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_w(x,y,z);
  return 0;
}

int MASA::masa_eval_p_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_p(x,y,z);
  return 0;
}

int MASA::masa_eval_rho_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_rho(x,y,z);
  return 0;
}

/* ------------------------------------------------
 *
 *         utility functions
 *
 * -----------------------------------------------
 */ 

int MASA::masa_get_name(string* name)
{
  masa_master_pointer->return_name(name); // set string to name
  return 0;
}

int MASA::masa_get_dimension(int* dim)
{
  masa_master_pointer->return_dim(dim); // set string to name
  return 0;
}

int MASA::masa_sanity_check()
{
  masa_master_pointer->sanity_check(); // set string to name
  return 0;
}

int MASA::masa_version_stdout()
{
  std::cout << "--------------------------------------------------------" << std::endl;
  std::cout << "MASA Library: Version = " << MASA_LIB_VERSION;
  std::cout << " (" << MASA::masa_get_numeric_version() << ")" << std::endl << std::endl;
  std::cout << "Build Date   = " << MASA_BUILD_DATE     << std::endl;
  std::cout << "Build Host   = " << MASA_BUILD_HOST     << std::endl;
  std::cout << "Build User   = " << MASA_BUILD_USER     << std::endl;
  std::cout << "Build Arch   = " << MASA_BUILD_ARCH     << std::endl;
  std::cout << "Build Rev    = " << MASA_BUILD_VERSION  << std::endl << std::endl;
  std::cout << "C++ Config   = " << MASA_CXX << " "     << MASA_CXXFLAGS << std::endl;
  //std::cout << "F90 Config   = " << MASA_FC MASA_FCFLAGS << std::endl;
  std::cout << "--------------------------------------------------------" << std::endl;
  return 0;
}

int MASA::masa_get_numeric_version()
{
  // Note: return format follows the versioning convention xx.yy.zz where
  //
  // xx = major version number
  // yy = minor version number
  // zz = micro version number
  //
  // For example:
  // v.   0.23  -> 002300 = 2300
  // v   0.23.1 -> 002301 = 2301
  // v. 10.23.2 -> 102302

  int major_version = 0;
  int minor_version = 0;
  int micro_version = 0;

  #ifdef MASA_MAJOR_VERSION
  major_version = MASA_MAJOR_VERSION;
  #endif

  #ifdef MASA_MINOR_VERSION
  minor_version = MASA_MINOR_VERSION;
  #endif

  #ifdef MASA_MICRO_VERSION
  micro_version = MASA_MICRO_VERSION;
  #endif

  return(major_version*10000 + minor_version*100 + micro_version);

}

