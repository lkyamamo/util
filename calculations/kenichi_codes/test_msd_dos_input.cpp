#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include <cstdlib>
#include <cmath>
#include <cassert>

struct MDFrame
{
    int natoms;
    float lattice[6];

    std::vector<float> x, y, z;
    std::vector<float> vx, vy, vz;
    std::vector<std::string> name;

	void print()
	{
		std::cout << " # of Atoms : " << natoms << std::endl;
		std::cout << " Lattice Consts : " <<
		lattice[0] << " " << lattice[1] << " " << lattice[2] << " " <<  
		lattice[3] << " " << lattice[4] << " " << lattice[5] << std::endl;

		for(int i=0; i<natoms; i++)
		{
			std::cout << name[i] << " " << x[i] << " " << y[i] << " " << z[i] << " " << 
			vx[i] << " " << vy[i] << " " << vz[i] << std::endl;
		}
	}
};

void print_mdframe(MDFrame &frame)
{
    std::cout << " # of Atoms : " << frame.natoms << std::endl;
    std::cout << " Lattice Consts : " <<
              frame.lattice[0] << " " << frame.lattice[1] << " " << frame.lattice[2] << "\n" <<
              frame.lattice[3] << " " << frame.lattice[4] << " " << frame.lattice[5];

    /*
    	for(int i=0; i<frame.natoms; i++)
    	{
    		std::cout << frame.name[i] << " " << frame.x[i] << " "
    							<< frame.y[i] << " " << frame.z[i] << std::endl;
    	}
    */
}

MDFrame read_single_mdframe(std::ifstream &in)
{
    MDFrame mdframe;
    std::string str;

    std::getline(in,str);
    mdframe.natoms = std::atoi(str.c_str());

    double dummy;
    std::stringstream ss;
    std::getline(in,str);
    ss << str;
    /*
    	ss >> mdframe.lattice[0] >> dummy >> dummy >>
    				dummy >> mdframe.lattice[1] >> dummy >>
    				dummy >> dummy >> mdframe.lattice[2];
    */
    ss >> mdframe.lattice[0] >> mdframe.lattice[1] >> mdframe.lattice[2] >>
	    mdframe.lattice[3] >> mdframe.lattice[4] >> mdframe.lattice[5];

    std::cout << mdframe.lattice[0] << mdframe.lattice[1] << mdframe.lattice[2] <<
    mdframe.lattice[3] << mdframe.lattice[4] << mdframe.lattice[5];

    mdframe.name.resize(mdframe.natoms);
    mdframe.x.resize(mdframe.natoms);
    mdframe.y.resize(mdframe.natoms);
    mdframe.z.resize(mdframe.natoms);
    mdframe.vx.resize(mdframe.natoms);
    mdframe.vy.resize(mdframe.natoms);
    mdframe.vz.resize(mdframe.natoms);

    for (int i=0; i<mdframe.natoms; i++)
    {
        std::string name;
	int id, idummy; 
        float x,y,z,vx,vy,vz, fdummy;

        std::stringstream ss;
        std::getline(in,str);
        ss << str;
        //ss >> name >> x >> y >> z >> vx >> vy >> vz;
        //ss >> name >> x >> y >> z >> id >> idummy >> vx >> vy >> vz; 
        //ss >> name >> x >> y >> z >> fdummy >> id >> vx >> vy >> vz; 
	//
#ifdef NNQMD
        ss >> name >> x >> y >> z >> dummy >> id;  
#else
        //ss >> name >> x >> y >> z >> dummy >> id >> vx >> vy >> vz;
        ss >> name >> x >> y >> z >> vx >> vy >> vz;
	//std::cout << name << " " << x << " " << y << " " << z << " " << vx << " " << vy << " " << vz << std::endl;
	id = i+1; 
#endif

        mdframe.name[id-1] = name;
        mdframe.x[id-1] = x;
        mdframe.y[id-1] = y;
        mdframe.z[id-1] = z;
        mdframe.vx[id-1] = vx;
        mdframe.vy[id-1] = vy;
        mdframe.vz[id-1] = vz;
        //std::cout << i << " " << id-1 << " " << name << " " << x << " " << y << " " << z << " " << vx << " " << vy << " " << vz << std::endl;
/*
        mdframe.name.push_back(name);
        mdframe.x.push_back(x);
        mdframe.y.push_back(y);
        mdframe.z.push_back(z);
        mdframe.vx.push_back(vx);
        mdframe.vy.push_back(vy);
        mdframe.vz.push_back(vz);
*/
    }

    return mdframe;
};