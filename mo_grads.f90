module class_GrADS
  implicit none
  private
  public :: GrADS, circle_area, circle_print

  real :: pi = 3.1415926535897931d0 ! Class-wide private constant

  type GrADS
  
	  INTEGER:: nx;
	  REAL(dp):: lonw;
	  REAL(dp):: dlon;
	  INTEGER:: ny;
	  REAL(dp):: lats;
	  REAL(dp):: dlat;
	  INTEGER:: nz;
	  string zmapping; // LINEAR or LEVELS
	  vector <REAL(dp)::> pz;		// z levels 
	  REAL(dp):: z;
	  REAL(dp):: dz;
	  INTEGER:: period;
	  string startDate; 	// e.g. 00:00Z20oct1997;
	  string interval;
	  INTEGER:: var_start;			// starting variable
	  INTEGER:: nvar;			// number of variable
	  vector <var> pvar;
	  vector <REAL::> data;	// binary data buffer
	  vector <var> pdata;		// data buffer variables  
     real :: radius
  end type GrADS
contains
  function circle_area(this) result(area)
    type(GrADS), intent(in) :: this
    real :: area
    area = pi * this%radius**2
  end function circle_area

  subroutine circle_print(this)
    type(GrADS), intent(in) :: this
    real :: area
    area = circle_area(this)  ! Call the circle_area function
    print *, 'GrADS: r = ', this%radius, ' area = ', area
  end subroutine circle_print
end module class_GrADS


#include <iostream>
using std:: endl;
using std:: cerr;
#include <fstream>
using std:: ofstream;
using std:: ifstream;
using std:: fstream;
using std:: ios;
#include <string>
using std::string;
#include <ctype.h>

#ifndef UNIX
#include <stdlib.h>
#define min __min
#else
#include <algorithm>
using std::min;
#endif
#ifdef UNIX
#include <strstream>
using std::istrstream;
#define istringstream istrstream
#else
#include <sstream>
using std::istringstream;
#endif

#ifndef UNIX
#define strcasecmp _stricmp
#endif

#include <vector>
using namespace std;

#include "math.h"

#include "constant.h"
#include "gradsclass.h"
#include "endianclass.h"

#ifndef GRADS
#define GRADS

#ifndef ofstream
#include <fstream>
// using namespace std;
using std:: ofstream;
using std:: ifstream;
using std:: fstream;
using std:: ios;
#endif

#ifndef string
#include <string>
using std::string;
#endif

#include <vector>
using namespace std;

#include "constant.h"
#ifndef _____ENDIAN_____
#include "endianclass.h"
#endif

class var {
public:
	string abrev;	// abrev[8]
//	char abrev[8];
	INTEGER:: levs;
	INTEGER:: units;
	//	char description[40];
	string description;	// description[40]
	REAL(dp):: a;
	REAL(dp):: b;
	INTEGER:: c;  // index for determining unit converion whether including StdV or not
	REAL(dp):: unitout(REAL(dp):: StdV, REAL(dp):: val);		// convert unit from internal to output, out=a*val+b
};


module class_GrADS
  implicit none
  private
  public :: GrADS, bin_file, ctldir, 

  real :: pi = 3.1415926535897931d0 ! Class-wide private constant

  type GrADS
class grads {
public:
	fstream bin_file;
	string ctldir;
	string fname;
	string map_fname;
	string title;
	REAL:: undef;
	string options;
	LOGICAL:: little_endian;	// true = data in little_endian format (default); false = data in big_endian format)
	LOGICAL:: yrev;
	LOGICAL:: zrev;
	LOGICAL:: ltemplate;	// template
	LOGICAL:: latlon;		// true= data in latlon system; false= utm system (default)
	INTEGER:: nx;
	REAL(dp):: lonw;
	REAL(dp):: dlon;
	INTEGER:: ny;
	REAL(dp):: lats;
	REAL(dp):: dlat;
	INTEGER:: nz;
	string zmapping; // LINEAR or LEVELS
	vector <REAL(dp)::> pz;		// z levels 
	REAL(dp):: z;
	REAL(dp):: dz;
	INTEGER:: period;
	string startDate; 	// e.g. 00:00Z20oct1997;
	string interval;
	INTEGER:: var_start;			// starting variable
	INTEGER:: nvar;			// number of variable
	vector <var> pvar;
	vector <REAL::> data;	// binary data buffer
	vector <var> pdata;		// data buffer variables
	grads (string fname0="", string title0="", REAL:: undef0=MISSING,
		INTEGER:: nx0=1, REAL(dp):: lon0=0, REAL(dp):: dlon0=1,
		REAL(dp):: utme0=0, REAL(dp):: dutme0=1,
		INTEGER:: ny0=1, REAL(dp):: lat0=0, REAL(dp):: dlat0=1,
		REAL(dp):: utmn0=0, REAL(dp):: dutmn0=1,
		INTEGER:: nz0=1, REAL(dp):: z0=0, REAL(dp):: dz0=100,
		INTEGER:: period0=0, string startDate0="", string interval0="1hr",
		INTEGER:: var_start0=0, INTEGER:: nvar0=0);
	~grads();
	// 	void assign(vector<var> pvar0);
	LOGICAL:: close();
    LOGICAL:: getOffset(REAL(dp):: x, REAL(dp):: y, INTEGER:: * offset);
	INTEGER:: imatrix(string varname );
	INTEGER:: ivar(INTEGER:: nmatrix );
	LOGICAL:: open();
	LOGICAL:: open_for_append();
	LOGICAL:: open_for_write();
//	LOGICAL:: read(REAL(dp):: x, REAL(dp):: y, string varname, REAL(dp):: * value,
//		REAL(dp):: minValue, REAL(dp):: maxValue, REAL(dp):: defaultValue);
	LOGICAL:: read(REAL(dp):: x, REAL(dp):: y, string varname, REAL(dp):: * value,
		REAL(dp):: minValue=FLT_SMALL, REAL(dp):: maxValue=FLT_HUGE, REAL(dp):: defaultValue=0.);
	LOGICAL:: read_ctl(string ctl_fname);
	LOGICAL:: readBuffer (REAL(dp):: x, REAL(dp):: y, INTEGER:: k, REAL(dp):: * value);
    LOGICAL:: readBinfile(string varname);
	LOGICAL:: readBinfile (REAL(dp):: x, REAL(dp):: y, string var, REAL(dp):: * value);
	void write_grid_ctl(string fname);
	//	void write_grid_ctl_2(string ctl_fname);
	void write_station_ctl(string ctl_fname);
};

class rpthdr {
public:
	char id[8]; 		  /* Character station id			*/
	REAL:: lat;			  /* Latitude of report 			*/
	REAL:: lon;			  /* Longitude of report			*/
	REAL:: t;			  /* Time in relative grid units	*/
	INTEGER::  nlev;			  /* Number of levels following 	*/
	INTEGER:: flag;			  /* Level independent var set flag */
	rpthdr();
	~rpthdr();
};
extern void UTMtoLonLat(REAL(dp)::, REAL(dp)::, REAL(dp):: *, REAL(dp):: *);
void do_tolower( string *words );
void do_toupper( string *words );
#endif


void do_tolower( string *words )
{
	
	string caps( "ABCDEFGHIJKLMNOPQRSTUVWXYZ" );
	
	string::size_type pos = 0;
	while (( pos = (*words).find_first_of( caps, pos )) != string::npos )
		(*words)[ pos ] = tolower( (*words)[pos] );
}
void do_toupper( string *words )
{
	string lowercase( "abcdefghijklmnopqrstuvwxyz" );
	string::size_type pos = 0;
	while (( pos = (*words).find_first_of( lowercase, pos )) != string::npos )
		(*words)[ pos ] = toupper( (*words)[pos] );
}

REAL(dp):: var:: unitout(REAL(dp):: StdV, REAL(dp):: val) {
	if ( 0 == c) return (a*val+b);
	else if ( 1==c) return (a*StdV*val+b);
	else return (val);
}

grads::  grads(string fname0, string title0, REAL:: undef0,
			   INTEGER:: nx0, REAL(dp):: lon0, REAL(dp):: dlon0,
			   REAL(dp):: utme0, REAL(dp):: dutme0,
			   INTEGER:: ny0, REAL(dp):: lat0, REAL(dp):: dlat0,
			   REAL(dp):: utmn0, REAL(dp):: dutmn0,
			   INTEGER:: nz0, REAL(dp):: z0, REAL(dp):: dz0,
			   INTEGER:: period0, string startDate0, string interval0,
			   INTEGER:: var_start0, INTEGER:: nvar0) {
	ctldir="";
	fname=fname0;
	map_fname=fname0+".map";
	title=title0.c_str();
	undef=undef0;
	options="";
	little_endian=true;
	yrev=false;
	zrev=false;
	ltemplate=false;
	latlon=false;
	nx=nx0;
	lonw=lon0;
	dlon=dlon0;
	utme=utme0;
	dutme=dutme0;
	ny=ny0;
	lats=lat0;
	dlat=dlat0;
	utmn=utmn0;
	dutmn=dutmn0;
	nz=nz0;
	zmapping="LINEAR";
	z=z0;
	dz=dz0;
	period=period0;
	startDate=startDate0;
	interval=interval0;
	var_start=var_start0;
	nvar=nvar0;
}

grads:: ~grads() {
}
/*
void grads:: assign(vector<var> pvar0) {
pvar.assign(pvar0.begin(),pvar0.end());
nvar=pvar.size();
}
*/
LOGICAL:: grads:: close(){
	// close binary file
	bin_file.close();
	return(true);
}
LOGICAL:: grads:: getOffset(REAL(dp):: x, REAL(dp):: y, INTEGER:: * offset){
	if (latlon) {	// latlon coordinate system
		REAL(dp):: lon,lat;
		UTMtoLonLat(x,y,&lon,&lat);
		// Check whether the range of ctl file cover the domain of simulation?
		if ((lonw-dlon/2. > lon)||
			(lonw+dlon*((REAL(dp)::)nx+0.5) < lon) ||
			(lats-dlat/2. > lat)||
			(lats+dlat*((REAL(dp)::)ny+0.5) < lat)	) {
#if DEBUG
			cerr<<"Domain of simulation is outside the range of the ctl file!"<<endl
				<<"Default value for those ranges."<<endl;
#endif
			return (false);
		}
		if ( yrev ) {
			*offset =
				+((INTEGER::)((lats+dlat*(ny+0.5)-lat)/dlat))*nx	//lat offset
				+((INTEGER::)((lon-(lonw-dlon/2.))/dlon)) +1				// lon offset
				-1;  // start from zero 
		}
		else {
			*offset =
				+((INTEGER::)((lat-(lats-dlat/2.))/dlat))*nx				//lat offset
				+((INTEGER::)((lon-(lonw-dlon/2.))/dlon))+1				//lon offset
				-1;  // start from zero  
		}
	}
	else {	// utm coordinate system
		// Check whether the range of ctl file cover the domain of simulation?
		if ((utme-dutme/2. > x)||
			(utme+dutme*((REAL(dp)::)nx+0.5) < x) ||
			(utmn-dutmn/2. > y)||
			(utmn+dutmn*((REAL(dp)::)ny+0.5) < y)	) {
#if DEBUG
			cerr<<"Domain of simulation is outside the range of the ctl file!"<<endl
				<<"Default value for those ranges is used."<<endl;
#endif
			return (false);
		}
		if ( yrev ) {
			*offset =
				+((INTEGER::)((utmn+dutmn*(ny+0.5)-y)/dutmn))*nx	//y offset
				+((INTEGER::)((x-(utme-dutme/2.))/dutme)) +1				// x offset
				-1;  // start from zero 
		}
		else {
			*offset =
				+( (INTEGER::)((y-(utmn-dutmn/2.))/dutmn) )*nx				//y offset
				+( (INTEGER::)((x-(utme-dutme/2.))/dutme) )+1 			//x offset
				-1;  // start from zero  
		}
	}
	return (true);
}
INTEGER:: grads:: imatrix(string varname ) {
	// to get the sequence of the matrix that the variable belong to  "varname"
	INTEGER:: imatrix =0;
	vector <var>::iterator it;
	for (it = pvar.begin(); it != pvar.end(); it++) {
		if (0==varname.compare(it->abrev)) {
			return (imatrix);
		}
		else {
			if (0 != it->levs )  { imatrix +=  it->levs; }
			else { imatrix ++; }
		}
	}
	cerr << "variable: "<<varname<<" does not exist in "<< fname<<endl;
	return(iMISSING);
}
INTEGER:: grads:: ivar(INTEGER:: nmatrix ) {
	// to get the sequence of the matrix that the variable belong to  "varname"
	INTEGER:: imatrix=0;
	INTEGER:: ivar =0;
	vector <var>::iterator it;
	for (ivar=0, it = pvar.begin(); it != pvar.end(); it++, ivar++) {
		if (0 != it->levs )  { imatrix +=  it->levs; }
		else { imatrix ++; }
		if (imatrix>nmatrix) return(ivar);
	}
	cerr << "data # "<<nmatrix<<" is out of the range."<<endl;
	return(iMISSING);
}

LOGICAL:: grads:: open(){
	// open binary file for read
	string fullname=ctldir+fname;
	bin_file.open(fullname.c_str(),ios::binary|ios::in);
	if (!bin_file) {
		cerr<<"Open "<<fullname.c_str()<<" error!"<<endl;
		exit(-1);
	}
	return (true);
}
LOGICAL:: grads:: open_for_append() {
	string fullname=ctldir+fname;
	bin_file.open(fullname.c_str(),ios::binary|ios::out|ios::app);
	if (!bin_file) {
		cerr<<"Open "<<fullname.c_str()<<" error!"<<endl;
		return(false);
	}
	return (true);
}
LOGICAL:: grads:: open_for_write() {
	string fullname=ctldir+fname;
	bin_file.open(fullname.c_str(),ios::binary|ios::out|ios::trunc);
	if (!bin_file) {
		cerr<<"Open "<<fullname.c_str()<<" error!"<<endl;
		return(false);
	}
	return (true);
}
LOGICAL:: grads:: read(REAL(dp):: x, REAL(dp):: y, string varname, REAL(dp):: * value,
	REAL(dp):: minValue, REAL(dp):: maxValue, REAL(dp):: defaultValue){
	REAL(dp):: tmp;
	LOGICAL:: in_data_buffer = false;
	LOGICAL:: ok=false;
	INTEGER:: k=0;	// maxtrix number
	vector <var>:: iterator it;
	for (it=pdata.begin();it!=pdata.end();it++) {
		if( 0==varname.compare(it->abrev) ) {
			in_data_buffer = true;
			break;
		}
		k++;
	}
	if ( in_data_buffer) {
		ok=readBuffer (x, y, k, &tmp);
	}
	else {
		readBinfile(varname);
		ok=readBuffer (x, y, k, &tmp);
	}
	if ((tmp > minValue) && (tmp < maxValue)) {
		*value=tmp;
	}
	else {
		*value=defaultValue; 
	}
	return (ok);
}
LOGICAL:: grads:: readBinfile(string varname){
	// read entire matrix of binary file of "varname" into data buffer
	extern Endian myComputer;
	LOGICAL:: read_ok = false;
	INTEGER:: errcount =0;
	INTEGER:: offset,k;
	REAL:: * pFloatValue;
	pFloatValue = new REAL::[nx*ny];
	for (INTEGER:: i = 0; i<nx*ny; i++) {
		pFloatValue[i]=(REAL::) MISSING;
	}
	if(iMISSING == (k=imatrix(varname)) ) {
		cerr<<varname<<" does not exist in "<<fname<<"!"<<endl;
		return (false);
	};
	while ((!read_ok) && (10>errcount)) {
		offset = nx*ny*k;	// get offset from file orig.
		bin_file.seekg(offset*sizeof(REAL::),ios::beg);
		bin_file.read((char *) pFloatValue,nx*ny*sizeof(REAL::));
		if (bin_file.good()) {
			for (INTEGER:: i = 0; i< nx*ny; i++) {
				myComputer.CheckAndSwapEndian((REAL:: *) &(pFloatValue[i]), little_endian);
				data.push_back(pFloatValue[i]);
			}
			vector <var>::iterator it;
			for (it = pvar.begin(); it != pvar.end(); it++) {
				if (0==varname.compare(it->abrev)) {
					pdata.push_back(*it);
					break;
				}
			}
			return (true);
		}
		else {
			errcount ++;
			cerr <<"read binfile error! try "<<errcount<<" times!"<< endl;
			if (bin_file.bad()) {
				cerr<<"critical error!"<<endl;
				exit(99);
			}
			else if(bin_file.fail()) {
				continue;
			}
		}		
	}
	delete [] pFloatValue;
	return (true);
}
LOGICAL:: grads:: readBinfile(REAL(dp):: x, REAL(dp):: y, string varname, REAL(dp):: * value){
	extern Endian myComputer;
	LOGICAL:: read_ok = false;
	INTEGER:: errcount =0;
	// read binary file to obtain the value of "varname" at (utme=x,utmn=y)
	INTEGER:: offset,k;
	REAL:: floatvalue=(REAL::) MISSING;
	if(iMISSING == (k=imatrix(varname)) ) {
		return (false);
	};
	// get offset from the variable orig.
	if (!getOffset(x,y,&offset) )  {
		return (false);
	}
	while ((!read_ok) && (10>errcount)) {
		offset += nx*ny*k;	 // get offset from file orig.
		bin_file.seekg(offset*sizeof(REAL::),ios::beg);
		bin_file.read((char *)&floatvalue,1*sizeof(REAL::));
		if (bin_file.good()) {
			myComputer.CheckAndSwapEndian((REAL:: *)&floatvalue, little_endian);
			if (undef == floatvalue) {return (false); }
			else {
				*value = (REAL(dp)::) floatvalue;
				read_ok=true;
				return (true);
			}
		}
		else {
			errcount ++;
			cerr <<"read binfile error! try "<<errcount<<" times!"<< endl;
			if (bin_file.bad()) {
				cerr<<"critical error!"<<endl;
				exit(99);
			}
			else if(bin_file.bad()) {
				continue;
			}
		}
	}
	return (true);
}

LOGICAL:: grads:: readBuffer(REAL(dp):: x, REAL(dp):: y, INTEGER:: k, REAL(dp):: * value){
	// read data buffer to obtain the value of "varname" at (utme=x,utmn=y)
	// k: number of matrix
	INTEGER:: offset;
	REAL:: * floatvalue= new REAL::;
	*floatvalue=(REAL::) MISSING;
	if (getOffset(x,y,&offset) )	// get offset from the variable orig.
	{
		floatvalue = data.begin()+k*nx*ny+offset;
		*value=(REAL(dp)::) *floatvalue;
		if (undef ==  *floatvalue) {
			*value = MISSING;
			return (false);
		}
		return (true);
	}
	else return (false);
}
LOGICAL:: grads:: read_ctl(string ctl_fname) {
	// Read description file (Ctl file) for grid data
	// return (false) for outside the range or missing value
	// otherwise return (true)
	const INTEGER:: lineSize=1024;
	char textline[lineSize];
	char tmpBuffer[40];
	REAL(dp):: insertz;
	var insert;
	string buffer, buffer2;
	string seperators("\\/");
	string:: size_type pos;

	ifstream f_ctl(ctl_fname.c_str(),ios::in);
	if (!f_ctl) {
		cerr<<"oops! unable to open file."
			<<ctl_fname<<"bailing out!\n";
		return(false);
	}
	while(f_ctl.getline(textline,lineSize,'\n')) {
		istringstream input_istring(textline,lineSize);
		input_istring>>buffer;
		do_toupper(&buffer);
		if ( 0 == buffer.compare("DSET"))  {
			input_istring>>fname;
			if ("^"==fname.substr(0,1)) {
				fname.erase(0,1);
			}
			else if ( (pos =fname.find_last_of(seperators) )
				!= string::npos ) {
				ctldir=fname.substr(0,pos+1);
				fname.erase(0,pos+1);
			}		
		}
		else if(0==buffer.compare("TITLE")) {
			input_istring.get(tmpBuffer,40);
			title=tmpBuffer;
		}
		else if (0==buffer.compare("UNDEF")) {
			input_istring>>undef;
		}
		else if (0==buffer.compare("OPTIONS")) {
			buffer2=textline;
			do_toupper(&buffer2);
			if ( (pos = buffer2.find("BIG_ENDIAN")) != string::npos) little_endian=false;
			if ( (pos = buffer2.find("LITTLE_ENDIAN")) != string::npos) little_endian=true;
			if ( (pos = buffer2.find("YREV")) != string::npos) yrev=true;
			if ( (pos = buffer2.find("ZREV")) != string::npos) zrev=true;
		}
		else if (0==buffer.compare("#OPTIONS")) {
			buffer2=textline;
			do_toupper(&buffer2);
			if ( (pos = buffer2.find("LATLON")) != string::npos) latlon=true;
			if ( (pos = buffer2.find("UTM")) != string::npos) latlon=false;;
		}
		else if (0==buffer.compare("XDEF")) {
			input_istring >> nx  >> buffer;
			do_toupper(&buffer);
			if (0==buffer.compare("LINEAR")) {
				input_istring>> lonw >> dlon;
			}
		}
		else if (0==buffer.compare("#XDEF")) {
			input_istring >> nx  >> buffer;
			do_toupper(&buffer);
			if (0==strcasecmp("LINEAR",buffer.c_str())) {
				input_istring>> utme >> dutme;
			}
		}
		else if (0==buffer.compare("YDEF")) {
			input_istring >> ny  >> buffer;
			do_toupper(&buffer);
			if ("LINEAR",buffer.c_str()) {
				input_istring>> lats >> dlat;
			}
		}
		else if (0==buffer.compare("#YDEF")) {
			input_istring >> ny  >> buffer;
			do_toupper(&buffer);
			if (0==buffer.compare("LINEAR")) {
				input_istring>> utmn >> dutmn;
			}
		}
		else if (0==buffer.compare("ZDEF")) {
			input_istring >> nz  >> buffer;
			do_toupper(&buffer);
			if (0==buffer.compare("LINEAR")) {
				zmapping="LINEAR";
				input_istring>> z >> dz;
			}
			else if (0==buffer.compare("LEVELS")) {
				zmapping="LEVELS";
				for (INTEGER:: i = 0; i< nz; i++) {
					input_istring>> insertz;
					pz.push_back(insertz);
				}
			}
		}
		else if (0==buffer.compare("TDEF")) {
			input_istring >> period  >> buffer;
			do_toupper(&buffer);
			if (0==buffer.compare("LINEAR")) {
				input_istring>> startDate >> interval;
			}
		}
		else if (0==buffer.compare("VARS")) {
			input_istring>>nvar;
			for (INTEGER:: i=0;i<nvar;i++) {
				f_ctl.getline(textline,lineSize,'\n');
				istringstream input_istring2(textline);
				input_istring2>>insert.abrev>>insert.levs>>insert.units;
				input_istring2.get(tmpBuffer,40);
				insert.description=tmpBuffer;
				pvar.push_back(insert);
			}
		}
		else if (0==buffer.compare("ENDVARS")) {
		}
	}
	f_ctl.close();
	return(true);
}
/*
void grads:: write_grid_ctl_2(string ctl_fname) {
// Write description file (Ctl file) for grid data
ofstream f_ctl; 
f_ctl.open(ctl_fname.c_str(),ios::out|ios::trunc);
f_ctl<<"DSET "<<fname<<endl
<<"TITLE "<<title<<endl
<<"UNDEF "<<undef<<endl
<<"OPTIONS little_endian"<<endl
<<"XDEF "<<nx<<" LINEAR "<<lonw<<" "<<dlon<<endl
<<"#XDEF "<<nx<<" LINEAR "<<utme<<" "<<dutme<<endl
<<"YDEF "<<ny<<" LINEAR "<<lats<<" "<<dlat<<endl
<<"#YDEF "<<ny<<" LINEAR "<<utmn<<" "<<dutmn<<endl
<<"ZDEF "<<nz<<" LINEAR "<<z<<" "<<dz<<endl
<<"TDEF "<<period<<" LINEAR "<<startDate<<" "<<interval<<endl
<<"VARS "<<nvar<<endl;
for (INTEGER:: i=0;i< nvar;i++){
f_ctl<<"v"<<i+1<<" 0 99 v"<<i+1<<endl; 
}
f_ctl<<"ENDVARS"<<endl;
f_ctl.close();
}
*/
void grads:: write_grid_ctl(string ctl_fname) {
	// Write description file (Ctl file) for grid data
	INTEGER:: nlevel;
	ofstream f_ctl;
	f_ctl.open(ctl_fname.c_str(),ios::out|ios::trunc);
	f_ctl<<"DSET ^"<<fname<<endl
		<<"TITLE "<<title<<endl
		<<"UNDEF "<<undef<<endl;
	if (little_endian) {
		options+=" little_endian";
	}
	else  {
		options+=" big_endian";
	}
	f_ctl<<"OPTIONS "<<options<<endl;
	if (latlon) {
		f_ctl<<"#OPTIONS latlon"<<endl;
	}
	else {
		f_ctl<<"#OPTIONS utm"<<endl;
	}
	f_ctl<<"XDEF "<<nx<<" LINEAR "<<lonw<<" "<<dlon<<endl
		<<"#XDEF "<<nx<<" LINEAR "<<utme<<" "<<dutme<<endl
		<<"YDEF "<<ny<<" LINEAR "<<lats<<" "<<dlat<<endl
		<<"#YDEF "<<ny<<" LINEAR "<<utmn<<" "<<dutmn<<endl;
	do_toupper(&zmapping);
	if (0==zmapping.compare("LINEAR")) {
		f_ctl<<"ZDEF "<<nz<<" LINEAR "<<z<<" "<<dz<<endl;
	}
	else if  (0==zmapping.compare("LEVELS")) {
		f_ctl<<"ZDEF "<<nz<<" LEVELS ";
		for (INTEGER:: i = 0; i<nz; i++) {
			f_ctl<<pz[i]<<" ";
		}
		f_ctl<<endl;
	}
	f_ctl<<"TDEF "<<period<<" LINEAR "<<startDate<<" "<<interval<<endl
		<<"VARS "<<nvar<<endl;
	for (INTEGER:: i=0;i<nvar;i++) {
		if (0==pvar[i].levs) {
			nlevel=0;
		}
		else  {
			nlevel = min(nz,pvar[i].levs);
		}
		f_ctl<<pvar[i].abrev.substr(0,8).c_str()<<"\t"<<nlevel<<" "<<pvar[i].units
			<<" "<<pvar[i].description<<endl;
	}
	f_ctl<<"ENDVARS"<<endl;
	f_ctl.close();
}
void grads:: write_station_ctl(string ctl_fname) {
	// Write description file (Ctl file) for station data
	INTEGER:: levs;
	ofstream f_ctl;
	f_ctl.open(ctl_fname.c_str(),ios::out|ios::trunc);
	f_ctl << "DSET ^"<<fname<<endl;
	f_ctl << "DTYPE station" <<endl;
	f_ctl << "STNMAP "<<map_fname<< endl;
	f_ctl << "TITLE "<<title<<endl;
	f_ctl << "UNDEF "<<undef<< endl;
	if (little_endian) {
		options+=" little_endian";
	}
	else  {
		options+=" big_endian";
	}
	f_ctl<<"OPTIONS "<<options<<endl;
	if (latlon) {
		f_ctl<<"#OPTIONS latlon"<<endl;
	}
	else {
		f_ctl<<"#OPTIONS utm"<<endl;
	}
	f_ctl << "TDEF "<<period<<" LINEAR "<<startDate<<" "<<interval<<endl;
	f_ctl << "VARS "<<nvar<<endl;
	for (INTEGER:: i=0;i<nvar;i++) {
		if (0==pvar[i].levs) {levs=0; } // surface variable
		else {levs=1;}		// level dependent variable
		f_ctl<<pvar[i].abrev.substr(0,8).c_str()<<"\t"<<levs<<" "<<pvar[i].units
			<<" "<<pvar[i].description<<endl;
	}
	f_ctl<<"ENDVARS"<<endl;
	f_ctl.close();
}
rpthdr:: rpthdr() {
	strcpy(id,"none");
	lat =0.;
	lon =0.;
	t = 0.;
	nlev=0;
	flag=1;
}
rpthdr:: ~rpthdr() {}

void UTMtoLonLat(REAL(dp):: utme, REAL(dp):: utmn, REAL(dp):: * lon, REAL(dp):: * lat)
{
	// convert UTM coordinate (m) to longitude/latitude (degree)
	//	REAL(dp):: *lat,*lon;		//	 latitude  &   longitude
	REAL(dp):: f0,y0,f1,f;
	REAL(dp):: a,b; 				// long radius	&  short radius
	REAL(dp):: v;					//	 m / sec  for latitude
	REAL(dp):: u;					//	 m / sec  for longitude
	REAL(dp):: delta_sec,delta_nata;
	REAL(dp):: r1,n1;				//	meridian radius  &	prime vertical radius
	REAL(dp):: esq; 				//	first flattening
	//	utme=275300.0;
	//	utmn=2765850.0;
	utme=utme-250000;
	f0=23.5/180.*PI;
	a=6378160.;
	b=6356774.5161;
	esq=0.006694605;
	delta_sec=sin(1.0/3600.0/180.0*PI);
	r1=a*(1-esq)/pow(1-esq*pow(sin(f0),2),1.5);
	n1=a/pow(1-esq*pow(sin(f0),2),0.5);
	v=(111133.34686*(f0*180.0/PI+1.0/3600.0) -16039.10038*sin(2.0*f0+2.0/3600.0/180.0*PI)+16.83394*sin(4.0*f0+4.0/3600.0/180.0*PI)-0.021654*sin(6.0*f0+6.0/3600.0/180.0*PI)) - (111133.34686*(f0*180.0/PI)-16039.10038*sin(2.0*f0)+16.83394*sin(4.0*f0)-0.021654*sin(6.0*f0));
	y0=(utmn-2600000)/v/3600.0/180.0*PI;
	f1=f0 + y0;
	f= f0 + y0; //-  (pow(utme,2.0)*tan(f1))/(2*r1*n1*delta_sec) +	pow(utme,4.0)*tan(f1)*(1+3*pow(tan(f1),2))/(24*r1*n1*n1*n1*delta_sec);
	delta_nata=utme/cos(f1)/(n1*delta_sec) - pow(utme,3)*(1/cos(f1))*pow(tan(f1),2)/(3*n1*n1*n1*delta_sec);
	*lat=f*180.0/PI;
	*lon=121.0 + delta_nata/3600.0;
	u=n1*cos(*lat/180.0*PI)*sin(1/3600.0/180.0*PI);
	return;
}