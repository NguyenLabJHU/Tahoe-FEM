#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <boost/unordered_set.hpp>

char *combineString(char *cstr, const char *str, int num, int width);

class Contact {
public:
  Contact();
  Contact(int p1, int p2, 
	  double p1x, double p1y, double p1z,
	  double p2x, double p2y, double p2z,
	  double r1, double r2,
	  double pene,
	  double td,
	  double cr,
	  double r0,
	  double e0,
	  double nf, double tf,
	  double cx, double cy, double cz,
	  double nx, double ny, double nz,
	  double tx, double ty, double tz,
	  double vs, double is)
    :ptcl_1(p1), ptcl_2(p2), 
     point1_x(p1x), point1_y(p1y), point1_z(p1z), 
     point2_x(p2x), point2_y(p2y), point2_z(p2z), 
     radius_1(r1), radius_2(r2), 
     penetration(pene), 
     tangt_disp(td), 
     contact_radius(cr), 
     R0(r0), E0(e0), 
     normal_force(nf), tangt_force(tf), 
     contact_x(cx), contact_y(cy), contact_z(cz), 
     normal_x(nx), normal_y(ny), normal_z(nz), 
     tangt_x(tx), tangt_y(ty), tangt_z(tz), 
     vibra_t_step(vs), impact_t_step(is) {}

  bool operator==(const Contact &other) const {
    return ( (ptcl_2 == other.ptcl_1 && ptcl_1 == other.ptcl_2) || (ptcl_1 == other.ptcl_1 && ptcl_2 == other.ptcl_2) );
  }
  
  friend std::size_t hash_value(const Contact &c) {
    boost::hash<int> hasher;
    return hasher(c.ptcl_1 * c.ptcl_2);
  }

public:
  int ptcl_1;
  int ptcl_2;
  double point1_x;
  double point1_y;
  double point1_z;
  double point2_x;
  double point2_y;
  double point2_z;
  double radius_1;
  double radius_2;
  double penetration;
  double tangt_disp;
  double contact_radius;
  double R0;
  double E0;
  double normal_force;
  double tangt_force;
  double contact_x;
  double contact_y;
  double contact_z;
  double normal_x;
  double normal_y;
  double normal_z;
  double tangt_x;
  double tangt_y;
  double tangt_z;
  double vibra_t_step;
  double impact_t_step;
};

int main(int argc, char *argv[])
{
  if(argc < 3) {
    std::cout << std::endl 
	      << "-- merge particle contact info from multiple processes and remove redundance --" << std::endl
	      << "  Usage: mergeContact  file_prefix  start_snap  end_snap  number_of_processes" << std::endl
	      << "Example: mergeContact  dep_contact  80  100  384" << std::endl << std::endl;
    return -1;
  }

  const int OWID = 16; 
  const int OPREC = 6; 
  
  int startSnap = atoi(argv[2]);
  int endSnap   = atoi(argv[3]);
  int totalProc = atoi(argv[4]);

  boost::unordered_set<Contact> allContact; // all contacts, no redundance

  ///////////////////////////////////////////////////////////
  for (int snapLoop = startSnap; snapLoop <= endSnap; ++snapLoop) {
    allContact.clear();

    char ostr[50], argsuf[50];
    strcpy(argsuf, argv[1]); strcat(argsuf, "_");
    combineString(ostr, argsuf, snapLoop, 3);
    std::ofstream ofs(ostr);
    if(!ofs) { std::cout << "ofstream error." << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    ///////////////////////////////////////////////////////////
    for (int procLoop = 0; procLoop < totalProc; ++procLoop) {
      char istr[50];
      combineString(istr, argsuf, snapLoop, 3);
      char csuf[10];
      combineString(csuf, ".p", procLoop, 5);
      strcat(istr, csuf);
      //std::cerr << istr << std::endl;

      std::ifstream ifs(istr);
      if(!ifs) { std::cout << "ifstream error." << std::endl; exit(-1); }
      int totalContact;
      std::string s;
      int ptcl_1;
      int ptcl_2;
      double point1_x;
      double point1_y;
      double point1_z;
      double point2_x;
      double point2_y;
      double point2_z;
      double radius_1;
      double radius_2;
      double penetration;
      double tangt_disp;
      double contact_radius;
      double R0;
      double E0;
      double normal_force;
      double tangt_force;
      double contact_x;
      double contact_y;
      double contact_z;
      double normal_x;
      double normal_y;
      double normal_z;
      double tangt_x;
      double tangt_y;
      double tangt_z;
      double vibra_t_step;
      double impact_t_step;
      
      ifs >> totalContact 
	  >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s
	  >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s;
      //std::cerr << totalContact << std::endl;
      ///////////////////////////////////////////////////////////
      for (int contactLoop = 0; contactLoop < totalContact; ++contactLoop) {
	ifs >> ptcl_1
	    >> ptcl_2
	    >> point1_x
	    >> point1_y
	    >> point1_z
	    >> point2_x
	    >> point2_y
	    >> point2_z
	    >> radius_1
	    >> radius_2
	    >> penetration
	    >> tangt_disp
	    >> contact_radius
	    >> R0
	    >> E0
	    >> normal_force
	    >> tangt_force
	    >> contact_x
	    >> contact_y
	    >> contact_z
	    >> normal_x
	    >> normal_y
	    >> normal_z
	    >> tangt_x
	    >> tangt_y
	    >> tangt_z
	    >> vibra_t_step
	    >> impact_t_step;
	allContact.insert(Contact(ptcl_1,ptcl_2,point1_x,point1_y,point1_z,point2_x,point2_y,point2_z,
				  radius_1,radius_2,penetration,tangt_disp,contact_radius,R0,E0,
				  normal_force,tangt_force,
				  contact_x,contact_y,contact_z,
				  normal_x,normal_y,normal_z,
				  tangt_x,tangt_y,tangt_z,vibra_t_step,impact_t_step));
      } // contactLoop
      ifs.close();
      remove(istr);
      //std::cerr << allContact.size() << std::endl;
    } // procLoop
    ofs << std::setw(OWID) << allContact.size() << std::endl;
    ofs << std::setw(OWID) << "ptcl_1"
	<< std::setw(OWID) << "ptcl_2"
	<< std::setw(OWID) << "point1_x"
	<< std::setw(OWID) << "point1_y"
	<< std::setw(OWID) << "point1_z"
	<< std::setw(OWID) << "point2_x"
	<< std::setw(OWID) << "point2_y"
	<< std::setw(OWID) << "point2_z"
	<< std::setw(OWID) << "radius_1"
	<< std::setw(OWID) << "radius_2"
	<< std::setw(OWID) << "penetration"
	<< std::setw(OWID) << "tangt_disp"
	<< std::setw(OWID) << "contact_radius"
	<< std::setw(OWID) << "R0"
	<< std::setw(OWID) << "E0"
	<< std::setw(OWID) << "normal_force"
	<< std::setw(OWID) << "tangt_force"
	<< std::setw(OWID) << "contact_x"
	<< std::setw(OWID) << "contact_y"
	<< std::setw(OWID) << "contact_z"
	<< std::setw(OWID) << "normal_x"
	<< std::setw(OWID) << "normal_y"
	<< std::setw(OWID) << "normal_z"
	<< std::setw(OWID) << "tangt_x"
	<< std::setw(OWID) << "tangt_y"
	<< std::setw(OWID) << "tangt_z"
	<< std::setw(OWID) << "vibra_t_step"
	<< std::setw(OWID) << "impact_t_step"
	<< std::endl;
    
    boost::unordered_set<Contact>::const_iterator it;    
    for (it = allContact.begin(); it != allContact.end(); ++it)
      ofs << std::setw(OWID) << it->ptcl_1       
	  << std::setw(OWID) << it->ptcl_2       
	  << std::setw(OWID) << it->point1_x     
	  << std::setw(OWID) << it->point1_y     
	  << std::setw(OWID) << it->point1_z     
	  << std::setw(OWID) << it->point2_x     
	  << std::setw(OWID) << it->point2_y     
	  << std::setw(OWID) << it->point2_z     
	  << std::setw(OWID) << it->radius_1     
	  << std::setw(OWID) << it->radius_2     
	  << std::setw(OWID) << it->penetration  
	  << std::setw(OWID) << it->tangt_disp   
	  << std::setw(OWID) << it->contact_radius
	  << std::setw(OWID) << it->R0           
	  << std::setw(OWID) << it->E0           
	  << std::setw(OWID) << it->normal_force 
	  << std::setw(OWID) << it->tangt_force  
	  << std::setw(OWID) << it->contact_x    
	  << std::setw(OWID) << it->contact_y    
	  << std::setw(OWID) << it->contact_z    
	  << std::setw(OWID) << it->normal_x     
	  << std::setw(OWID) << it->normal_y     
	  << std::setw(OWID) << it->normal_z     
	  << std::setw(OWID) << it->tangt_x      
	  << std::setw(OWID) << it->tangt_y      
	  << std::setw(OWID) << it->tangt_z      
	  << std::setw(OWID) << it->vibra_t_step 
	  << std::setw(OWID) << it->impact_t_step
	  << std::endl;

    ofs.close();
  } // snapLoop

  return 0;
}

char *combineString(char *cstr, const char *str, int num, int width) {
  std::string obj(str);
  std::stringstream ss;
  ss << std::setw(width) << std::setfill('0') << std::right << num;
  obj += ss.str();
  return strcpy( cstr, obj.c_str() );
}


