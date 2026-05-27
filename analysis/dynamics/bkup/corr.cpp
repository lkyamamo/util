#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <set>
#include <regex>

#include <cstdlib>
#include <cmath>
#include <cassert>

#include <boost/histogram.hpp>
#include <boost/format.hpp>

using namespace boost::histogram;

namespace fs = std::filesystem;

struct bond_stat
{
	double d0d0 = 0.0;
	std::vector<double> d0dt; 

	std::vector<double> p1, p2; 
};

struct dist 
{
	double dx,dy,dz;
	dist(double _dx, double _dy, double _dz) : dx(_dx), dy(_dy), dz(_dz){};
};

typedef std::tuple<double, int, std::string, dist, bool, int, std::string> tuple_t;
typedef std::vector<std::vector<tuple_t>> nbrlist_t;

#define N2_CUTOFF 2.46
#define N3_CUTOFF 2.73

#define MAX_NEIGHBOR_MOLS 5
//#define HBOND_OUTER_CUTOFF 2.73
#define HBOND_OUTER_CUTOFF 4.0
#define HBOND_INNER_CUTOFF 1.3


bool nbrlist_compare(nbrlist_t &t1, nbrlist_t &t2)
{
	// return true if two nbrlists are empty 
	if(t1.size()==0 && t2.size()==0) return true;

	// check if number of atoms are the same
	if (t1.size() != t2.size()) return false;

	std::set<int> set1, set2;

	for(int i = 0; i<t1.size(); i++)
	{
		// check if number of nbrs of i-th atom are the same
		if(t1[i].size() != t2[i].size()) return false; 

		// check if the nbr global id is the same
		for(int j = 0; j < t1[i].size(); j++)
		{
			set1.insert(std::get<1>(t1[i][j]));
			set2.insert(std::get<1>(t2[i][j]));
		}
		if(set1 != set2) return false; 
	}

	return true; 
}

// to read .bnd file
typedef std::vector<std::vector<int>> simple_nbrlist_t;

simple_nbrlist_t get_simple_nbrlists(std::string const & filename)
{
	std::ifstream fin(filename);
	simple_nbrlist_t simple_nbrlists; 

	std::string buffer;

	int num_atoms = 0; 
	while(std::getline(fin, buffer)) {num_atoms++;};

	simple_nbrlists.resize(num_atoms);

	fin.clear();
	fin.seekg(0);

	std::stringstream ss;
	while(std::getline(fin, buffer))
	{
		ss << buffer;

		int id,iid,jid,type_,num_nbrs; 
		double x_,y_,z_,dr;
		ss >> id >> x_ >> y_ >> z_ >> type_ >> num_nbrs;

		iid = id-1;
		for(int i=0; i<num_nbrs; i++)
		{
			ss >> jid >> dr; 
			int jid1 = jid - 1; // zero-indexed
			simple_nbrlists[iid].push_back(jid1);
		}

		//std::cout << buffer << std::endl;; 
		ss.str(""); ss.clear();
	}

	return simple_nbrlists;
}

struct MDFrame
{
    std::string filename;

    int natoms;
    double lattice[6];

    std::vector<double> x, y, z;
    std::vector<double> vx, vy, vz;
    std::vector<int> mol_id;  // molecular Id
    std::vector<std::string> name;

    double apply_pbc(double x, double lattice)
    {
        if(x>=0.5*lattice) x-=lattice;
        if(x<-0.5*lattice) x+=lattice;
        return x;
    };

    void print()
    {
      std::cout << " # of Atoms : " << natoms << std::endl;
      std::cout << " Lattice Consts : " <<
              lattice[0] << " " << lattice[1] << " " << lattice[2] << " " <<
              lattice[3] << " " << lattice[4] << " " << lattice[5] << std::endl;
    }
};

struct NeighborList
{
	int max_neighbors;
	int num_atoms; 

	nbrlist_t nbrlist;
	nbrlist_t nbrlist_ncenter; 
	nbrlist_t nbrlist_hcenter; 

	nbrlist_t nbrlist_ncenter_long; 

	simple_nbrlist_t slist;

	NeighborList() {}; 

	NeighborList(MDFrame & mdframe, int _max=MAX_NEIGHBOR_MOLS) : max_neighbors(_max), num_atoms(mdframe.natoms)
	{
		std::string bndfile = std::regex_replace(mdframe.filename, std::regex("xyz"), "bnd");

		if(fs::exists(bndfile)) 
		{
			std::cout << mdframe.filename << " " << bndfile << " " << fs::exists(bndfile) << std::endl;
			slist = get_simple_nbrlists(bndfile); 
		} else {
			std::cout << mdframe.filename << " " << fs::exists(bndfile) << std::endl;
		}

		nbrlist.resize(mdframe.natoms);	
		nbrlist_ncenter.resize(mdframe.natoms);	
		nbrlist_hcenter.resize(mdframe.natoms);	

		nbrlist_ncenter_long.resize(mdframe.natoms);	


		for (int i=0; i<mdframe.natoms; i++)
		{
			//std::cout << mdframe.name[i] << " " << mdframe.x[i] << " " << mdframe.y[i] << " " << mdframe.z[i] << std::endl; 

			//// construct nbrlist from only O
			//if (mdframe.name[i].find("O") == std::string::npos) continue; 
			
			bool is_i_oxygen = mdframe.name[i].find("O") != std::string::npos;

			std::vector<int> idlist; 
			if(slist.size()==mdframe.natoms)
			{
				for(int j = 0; j<slist[i].size(); j++) idlist.push_back(slist[i][j]);
			} else {
				for(int i = 0; i<mdframe.natoms; i++) idlist.push_back(i);
			}

/*
			for (int j1=0; j1<idlist.size(); j1++)
				std::cout << idlist[j1] << " "; 
				std::cout << std::endl;
*/
			for (int j1=0; j1<idlist.size(); j1++)
			{
				int j = idlist[j1];
				if(i==j) continue;

				//std::cout << j << " ";

				bool is_j_oxygen = mdframe.name[j].find("O") != std::string::npos;

				// do not include WCs
				if (mdframe.name[j].find("X") != std::string::npos) continue;

				double dx = mdframe.x[j] - mdframe.x[i];
				double dy = mdframe.y[j] - mdframe.y[i];
				double dz = mdframe.z[j] - mdframe.z[i];

				dx = mdframe.apply_pbc(dx, mdframe.lattice[0]);
				dy = mdframe.apply_pbc(dy, mdframe.lattice[1]);
				dz = mdframe.apply_pbc(dz, mdframe.lattice[2]);

				double dr = sqrt(dx*dx + dy*dy + dz*dz); 
				if (dr == 0.0) continue;

				auto data = std::make_tuple(dr, j, mdframe.name[j], dist(dx,dy,dz), true, i, mdframe.name[i]);

				nbrlist[i].push_back(data);

				//// ignore if an atomic pair distance is beyond than h-bond outer cutoff
				//if (HBOND_OUTER_CUTOFF < dr) continue;

				if(dr < 1.2 && is_i_oxygen && not is_j_oxygen) nbrlist_ncenter[i].push_back(data);
				if(dr < 4.0 && not is_i_oxygen && is_j_oxygen) nbrlist_hcenter[i].push_back(data);
				if(dr < 2.46 && is_i_oxygen && not is_j_oxygen) nbrlist_ncenter_long[i].push_back(data);
			}
		}

		// sort hcenter list by distance
		for (std::vector<tuple_t> & n : nbrlist_hcenter)
			std::sort(begin(n), end(n), [](auto const & t1, auto const & t2) { return std::get<0>(t1) < std::get<0>(t2);} );
		// sort ncenter list by id
		for (std::vector<tuple_t> & n : nbrlist_ncenter)
			std::sort(begin(n), end(n), [](auto const & t1, auto const & t2) { return std::get<1>(t1) < std::get<1>(t2);} );
		// sort ncenter_long list by distance
		for (std::vector<tuple_t> & n : nbrlist_ncenter_long)
			std::sort(begin(n), end(n), [](auto const & t1, auto const & t2) { return std::get<0>(t1) < std::get<0>(t2);} );
	}

	void print()
	{
		for (int i=0; i<nbrlist_hcenter.size(); i++)
		{
			if (nbrlist_hcenter[i].size() == 0) continue; 
			std::cout << std::endl << i << "th neighbr list size : " << nbrlist_hcenter[i].size() << std::endl << std::endl;

			for (const auto & l : nbrlist_hcenter[i])
				std::cout << i << " " << std::get<0>(l) << " " << std::get<1>(l) << " " << std::get<2>(l) << std::endl; 
			std::cout << std::endl; 
		}
	}

};


MDFrame read_single_mdframe(std::ifstream &in, std::string _filename="NA")
{
    MDFrame mdframe;
    std::string str;

    std::getline(in,str);
    mdframe.filename = _filename;
    mdframe.natoms = std::atoi(str.c_str());

    double dummy;
    std::stringstream ss;
    std::getline(in,str);
    ss << str;

    ss >> mdframe.lattice[0] >> mdframe.lattice[1] >> mdframe.lattice[2];
    mdframe.lattice[3] = mdframe.lattice[4] = mdframe.lattice[5] = 90.0;

    mdframe.name.resize(mdframe.natoms);
    mdframe.x.resize(mdframe.natoms);
    mdframe.y.resize(mdframe.natoms);
    mdframe.z.resize(mdframe.natoms);
    mdframe.mol_id.resize(mdframe.natoms);
    //std::cout << mdframe.natoms << std::endl;

    for (int i=0; i<mdframe.natoms; i++)
    {
        std::string name;
        float x,y,z,vx,vy,vz;
        int id;

        std::stringstream ss;
        std::getline(in,str);
        ss << str;

#ifdef NNQMD
        ss >> name >> x >> y >> z >> dummy >> id;  
#else
        //ss >> name >> x >> y >> z >> dummy >> id >> vx >> vy >> vz;
        ss >> name >> x >> y >> z; 
	id = i+1; 
#endif
        mdframe.name[id-1] = name;
        mdframe.x[id-1] = x;
        mdframe.y[id-1] = y;
        mdframe.z[id-1] = z;
        mdframe.mol_id[id-1] = id;
    }

    return mdframe;
};

struct hbond_populations
{
	struct hb_population
	{
		int num_atoms, num_frames; 
		double C1, C2;
		//std::vector<bool **> hbond_maps;
		//std::vector<float **> hbond_maps_dr;
		std::vector<std::vector<std::vector<bool>>> hbond_maps;

		void init(double _c1, double _c2, int _nf, int _na)
		{
			C1 = _c1; C2 = _c2;
			num_frames = _nf; num_atoms = _na;
		};

		void allocate_new_hbond_map()
		{
			//bool ** hbond_map;
			std::vector<std::vector<bool>> hbond_map;
			//float ** hbond_map_dr;

			//hbond_map = new bool * [num_frames];
			hbond_map.resize(num_frames);
			//hbond_map_dr = new float * [num_frames];

			for(int i = 0; i < num_frames; i++)
			{
				//hbond_map[i] = new bool [num_atoms*num_atoms];
				hbond_map[i].resize(num_atoms*num_atoms);
				for(int j = 0; j < num_atoms*num_atoms;  j++) hbond_map[i][j] = false;

				//hbond_map_dr[i] = new float [num_atoms*num_atoms];
				//for(int j = 0; j < num_atoms*num_atoms;  j++) hbond_map_dr[i][j] = 0.0;
			}

			hbond_maps.push_back(hbond_map);
			//hbond_maps_dr.push_back(hbond_map_dr);
		}

		void sample(int tseries_id, int step, int jgid, int kgid, float dr)
		{
			//hbond_map[step][jgid*num_atoms + kgid] = true;
			hbond_maps[tseries_id][step][jgid*num_atoms + kgid] = true;
			//hbond_maps_dr[tseries_id][step][jgid*num_atoms + kgid] = dr;

			//std::cout << "sample: step,jgid,kgid:  " << step << " " << jgid << " " << kgid << std::endl;
		};

		std::vector<double> summary(void)
		{
			std::vector<double> corr(num_frames, 0.0);

			// get number of hbonded molecules at t = 0.
			for(int ts = 0; ts < hbond_maps.size(); ts++)
			{
				long int num_hb0 = 0;
				for(int idx = 0; idx < num_atoms*num_atoms; idx++)
					if(hbond_maps[ts][0][idx] == true) num_hb0++;
				//std::cout << "tseries_id, initial hbonds : " << ts << "/" << hbond_maps.size() << " " << num_hb0 << std::endl;

				corr[0] += num_hb0; // # of hbond in first step 
			}

			for(int ts = 0; ts < hbond_maps.size(); ts++)
			{
				for(int step = 1; step < num_frames; step++)
				{
					int counter = 0;
					//std::cout << "step/num_frames,counter " << step << " / " << num_frames << " " << counter << std::endl;
					for(int idx = 0; idx < num_atoms*num_atoms; idx++)
					{
						if(hbond_maps[ts][step-1][idx] == false) hbond_maps[ts][step][idx] = false;
						if(hbond_maps[ts][0][idx] == true && hbond_maps[ts][step][idx] == true) counter++;
					}
					corr[step] += counter; 
				}
			}

			std::cout << "corr.size() " << corr.size() << " done." << std::endl;

			return corr;
		}

	};

	hb_population N2_s, N3_s, N2_s2l; 
	hb_population N2_25, N3_25;
	hb_population N2_26, N3_26;

	void initialize(int num_frames, int num_atoms)
	{
		N2_s.init(0.0, N2_CUTOFF, num_frames, num_atoms);
		//N3_s.init(0.0, N3_CUTOFF, num_frames, num_atoms);
		//N2_s2l.init(N2_CUTOFF, N3_CUTOFF, num_frames, num_atoms);

	}

	void allocate_new_hbond_map_all()
	{
		N2_s.allocate_new_hbond_map();
		//N3_s.allocate_new_hbond_map();
		//N2_s2l.allocate_new_hbond_map();
		
		std::cout << "adding new hbmap " << get_num_tseries() << std::endl; 
	}

	int get_num_tseries()
	{
		//assert(N2_s.hbond_maps.size() == N3_s.hbond_maps.size());
		return N2_s.hbond_maps.size();
	}


};


struct hb_profiler
{
	std::vector<std::vector<std::vector<double>>> hb_profile;
	std::vector<std::vector<double>> hb_ave_counts;

	std::vector<std::vector<double>> dhb_ave; 

	int N_count_lt_24 = 0; 
	int N_count_lt_246 = 0; 
	int N_count_lt_28 = 0; 
	int N_count_samples = 0;
	
	void initialize(int num_atoms, std::vector<fs::path> filepaths)
	{
		int num_frames = filepaths.size();
		hb_profile.resize(num_frames);
		hb_ave_counts.resize(num_frames);

		for(int iframe = 0; iframe < num_frames; iframe++) 
		{
			hb_ave_counts[iframe].resize(2);

			hb_profile[iframe].resize(num_atoms);
			for (int iatom = 0; iatom < num_atoms; iatom++)
				hb_profile[iframe][iatom].resize(3);
		}
		dhb_ave.resize(num_frames);
	};

	void get_average_hb_distance(int iframe, NeighborList const & nbr)
	{
		for (int iatom = 0; iatom < nbr.nbrlist_ncenter_long.size(); iatom++)
		{
			const auto n = nbr.nbrlist_ncenter_long[iatom]; 

			if(n.size() >= 4) 
			{
				// H4 : hydrogen bonded H from N
				assert(std::get<2>(n[3]).find("H") != std::string::npos);

				const double dnh = std::get<0>(n[3]);
				if (1.8 < dnh && dnh < 2.46) dhb_ave[iframe].push_back(dnh);

				/*
				for(int ia=0; ia<4; ia++)
				{
					std::cout << std::get<0>(n[ia]) << " ";
					assert(std::get<2>(n[ia]).find("H") != std::string::npos);
				}
				std::cout << " dnh " << dnh << std::endl;
				*/
			}

		}
	}

	void get_neighbor_n_distances(int iframe, NeighborList const & nbr)
	{
		for (int iatom = 0; iatom< nbr.nbrlist_ncenter.size(); iatom++)
		{
			const auto n = nbr.nbrlist_ncenter[iatom]; 


			if(n.size() >= 3) 
			{
				int id_h1 = std::get<1>(n[0]);
				int id_h2 = std::get<1>(n[1]);
				int id_h3 = std::get<1>(n[2]);

				//std::cout << jgid << " " << kgid << " " << lgid << " " << 
				//std::get<2>(n[0]) << " " << std::get<2>(n[1]) << " " << std::get<2>(n[2]) << std::endl;

				assert(std::get<2>(n[0]).find("H") != std::string::npos);
				assert(std::get<2>(n[1]).find("H") != std::string::npos);
				assert(std::get<2>(n[2]).find("H") != std::string::npos);

				const auto nh1 = nbr.nbrlist_hcenter[id_h1]; 
				const auto nh2 = nbr.nbrlist_hcenter[id_h2]; 
				const auto nh3 = nbr.nbrlist_hcenter[id_h3]; 

				/*
				assert(nh1.size()>1);
				assert(nh2.size()>1);
				assert(nh3.size()>1);
				*/
				if(nh1.size()<1) continue;
				if(nh2.size()<1) continue;
				if(nh3.size()<1) continue;

				/*
				assert(std::get<2>(nh1[1]).find("N") != std::string::npos);
				assert(std::get<2>(nh2[1]).find("N") != std::string::npos);
				assert(std::get<2>(nh3[1]).find("N") != std::string::npos);
				*/
				if(std::get<2>(nh1[1]).find("O") == std::string::npos) continue;
				if(std::get<2>(nh2[1]).find("O") == std::string::npos) continue;
				if(std::get<2>(nh3[1]).find("O") == std::string::npos) continue;

				const double nj_norm = std::get<0>(nh1[1]);
				const double nk_norm = std::get<0>(nh2[1]);
				const double nl_norm = std::get<0>(nh3[1]);

				if(nj_norm < 2.46) N_count_lt_246++; 
				if(nk_norm < 2.46) N_count_lt_246++; 
				if(nl_norm < 2.46) N_count_lt_246++; 

				if(nj_norm < 2.8) N_count_lt_28++; 
				if(nk_norm < 2.8) N_count_lt_28++; 
				if(nl_norm < 2.8) N_count_lt_28++; 
				N_count_samples += 1;

				//std::cout << "N-H distance: " << iframe << " " << iatom << " " << nj_norm  << " " << nk_norm << " " << nl_norm  << std::endl;
				hb_profile[iframe][iatom] = {nj_norm, nk_norm, nl_norm};

				/*
				auto & hb = hb_profile[iframe][iatom];
				std::cout << "N-H distance: " << iframe << " " << iatom << " " << hb[0] << " " << hb[1] << " " << hb[2] 
					<< " " << hb_profile.size() << " " << hb_profile[iframe].size() << " " << hb_profile[iframe][iatom].size() << std::endl;
				*/

			}

		}

		hb_ave_counts[iframe] = {(float) N_count_lt_246/N_count_samples, (float) N_count_lt_28/N_count_samples};
		//std::cout << "N-H ave count: " << iframe << " " << hb_ave_counts[iframe][0] << " " << hb_ave_counts[iframe][1] << std::endl;
	}

	void save()
	{
		std::ofstream fout_hb_profile(std::string("hb_profile.dat"));	
		for (int iframe = 0; iframe < hb_profile.size(); iframe++)
		{
			//std::cout << iframe << std::endl;
			fout_hb_profile << iframe << " " << hb_ave_counts[iframe][0] << " " << hb_ave_counts[iframe][1] << std::endl;

			for (int iatom = 0; iatom < hb_profile[iframe].size(); iatom++)
			{
				auto & hbs = hb_profile[iframe][iatom];
				//std::cout << hbs[0] << " " << hbs[1] << " " << hbs[2] << " ";
				fout_hb_profile << hbs[0] << " " << hbs[1] << " " << hbs[2] << " ";
			}
			//std::cout << std::endl;
			fout_hb_profile << std::endl;
		}

		long count_total = 0;
		double mean_total = 0.0;

		for (int iframe = 0; iframe < hb_profile.size(); iframe++)
		{
			int count = dhb_ave[iframe].size();
			double mean = 0.0;
			double stdev = 0.0;

			if (count > 0)
			{
				for(int i = 0; i<count; i++) mean +=dhb_ave[iframe][i]; 
				mean /= count; 
				for(int i = 0; i<count; i++) stdev += std::pow((dhb_ave[iframe][i] - mean),2);
				stdev = std::sqrt(stdev/count);

				for(int i = 0; i<count; i++) mean_total +=dhb_ave[iframe][i]; 
				count_total += count; 
			}

			std::cout << "iframe,count,mean,stddev: " << iframe << " " << count << " " << mean << " " << stdev << std::endl;
		}

		mean_total /= count_total; 
		double stdev_total = 0.0;
		for (int iframe = 0; iframe < hb_profile.size(); iframe++)
			for(int i = 0; i < dhb_ave[iframe].size(); i++) stdev_total += std::pow((dhb_ave[iframe][i] - mean_total),2);

		stdev_total = std::sqrt(stdev_total/count_total);
		std::cout << "TOTAL: count,mean,stddev: " << count_total << " " << mean_total << " " << stdev_total << std::endl;

	}

};

struct rotational_correlations
{

	typedef std::map<int, std::vector<double>> corr_ref_t;
	std::vector<corr_ref_t> e1_refs, e2_refs, e3_refs; 

	int corr_length, corr_interval, total_tseries;
	float dt_per_frame; 

	struct corr_t
	{
		double value = 0.0; 
		long num_samples = 0; 

		void add(double v, long n=1)
		{
			value += v;
			num_samples +=n;
		};
	};

	std::vector<corr_t> e1_corr, e2_corr, e3_corr;
	std::vector<corr_t> e1_l2_corr, e2_l2_corr, e3_l2_corr;

	void initialize(int _corr_length, int _corr_interval, float _dt_per_frame, int _total_tseries)
	{
		corr_length = _corr_length; 
		corr_interval = _corr_interval; 
		dt_per_frame = _dt_per_frame; 
		total_tseries =  _total_tseries;

		e1_corr.resize(corr_length);
		e2_corr.resize(corr_length);
		e3_corr.resize(corr_length);

		e1_l2_corr.resize(corr_length);
		e2_l2_corr.resize(corr_length);
		e3_l2_corr.resize(corr_length);
	}

	void normalize_vector(std::vector<double> &v)
	{
		double dr_norm = 0.0;
		for (int i=0; i<3; i++) dr_norm += v[i]*v[i];
		dr_norm = std::sqrt(dr_norm);
		assert(dr_norm>0.0);
		for (int i=0; i<3; i++) v[i] /= dr_norm;
	}

	auto get_normalized_vectors(const std::vector<tuple_t> & n)
	{
		// H1 
		int jgid = std::get<1>(n[0]); 
		const auto nj_r = std::get<3>(n[0]);
		const double nj_norm = std::sqrt(nj_r.dx*nj_r.dx + nj_r.dy*nj_r.dy + nj_r.dz*nj_r.dz);
	
		// H2 
		int kgid = std::get<1>(n[1]);
		const auto nk_r = std::get<3>(n[1]);
		const double nk_norm = std::sqrt(nk_r.dx*nk_r.dx + nk_r.dy*nk_r.dy + nk_r.dz*nk_r.dz);

		// H3 
		int lgid = std::get<1>(n[2]);
		const auto nl_r = std::get<3>(n[2]);
		const double nl_norm = std::sqrt(nl_r.dx*nl_r.dx + nl_r.dy*nl_r.dy + nl_r.dz*nl_r.dz);

		std::vector<double> e1 = {nj_r.dx+nk_r.dx+nl_r.dx, nj_r.dy+nk_r.dy+nl_r.dy, nj_r.dz+nk_r.dz+nl_r.dz};
		std::vector<double> e2 = {nj_r.dx, nj_r.dy, nj_r.dz}; 
		std::vector<double> e3 = {nj_r.dx-nk_r.dx, nj_r.dy-nk_r.dy, nj_r.dz-nk_r.dz};

		normalize_vector(e1);
		normalize_vector(e2);
		normalize_vector(e3);

		return std::make_tuple(e1,e2,e3);
	}

	void register_ref_vectors(nbrlist_t const & nbrlist_ncenter)
	{
		corr_ref_t e1_ref, e2_ref, e3_ref;

		if(e1_refs.size() >= total_tseries) return;

		int num_atoms = nbrlist_ncenter.size();
		std::cout << "registering new reference in corr" << std::endl;

		for (int iatom = 0; iatom< nbrlist_ncenter.size(); iatom++)
		{
			///std::cout << "i,nbrlist_hcenter.size() " << iatom <<  " " << nbr.nbrlist_hcenter[iatom].size() << std::endl;
	
			auto const & n = nbrlist_ncenter[iatom]; 
	
			if(n.size() >= 3) 
			{
				auto vectors = get_normalized_vectors(n);
				e1_ref[iatom] = std::get<0>(vectors);
				e2_ref[iatom] = std::get<1>(vectors);
				e3_ref[iatom] = std::get<2>(vectors);
			}
		}

		e1_refs.push_back(e1_ref);
		e2_refs.push_back(e2_ref);
		e3_refs.push_back(e3_ref);

	};

	void compute_correlations(int iframe, nbrlist_t const & nbrlist_ncenter)
	{
		for (int ts = 0; ts < e1_refs.size(); ts++)
		{
			int local_time  = iframe - ts*corr_interval; 
			assert(local_time >= 0);
	
			if(local_time < corr_length) 
			{ 
				std::cout << "ts,local_step,iframe,clength,cinterval: " 
					<< ts << " " <<  local_time << " " <<  iframe << " "
					<< corr_length << " " << corr_interval << std::endl;

			} else {
				std::cout << "ts,iframe,clength,cinterval: " 
					<< ts << " finished "
					<< corr_length << " " << corr_interval << std::endl;
			}

			if(local_time >= corr_length) continue;

			corr_t e1_dtd0,  e2_dtd0, e3_dtd0;
	
			for (int ia = 0; ia< nbrlist_ncenter.size(); ia++)
			{
				auto const & n = nbrlist_ncenter[ia]; 
		
				if(n.size() >= 3) 
				{

					/*
					if (e1_refs[ts].find(ia) == e1_refs[ts].end()) std::cout << "WARIING: key " << ia << " doesn't exist in e1_refs" << std::endl; 
					if (e2_refs[ts].find(ia) == e2_refs[ts].end()) std::cout << "WARIING: key " << ia << " doesn't exist in e2_refs" << std::endl; 
					if (e3_refs[ts].find(ia) == e3_refs[ts].end()) std::cout << "WARIING: key " << ia << " doesn't exist in e3_refs" << std::endl; 
					*/
					bool e1_flag = e1_refs[ts].find(ia) != e1_refs[ts].end();
					bool e2_flag = e2_refs[ts].find(ia) != e2_refs[ts].end();
					bool e3_flag = e3_refs[ts].find(ia) != e3_refs[ts].end();

					if(e1_flag && e2_flag && e3_flag)
					{
						auto vectors = get_normalized_vectors(n);
						auto e1 = std::get<0>(vectors);
						auto e2 = std::get<1>(vectors);
						auto e3 = std::get<2>(vectors);

						e1_dtd0.add(e1_refs[ts][ia][0]*e1[0] + e1_refs[ts][ia][1]*e1[1] + e1_refs[ts][ia][2]*e1[2]);
						e2_dtd0.add(e2_refs[ts][ia][0]*e2[0] + e2_refs[ts][ia][1]*e2[1] + e2_refs[ts][ia][2]*e2[2]);
						e3_dtd0.add(e3_refs[ts][ia][0]*e3[0] + e3_refs[ts][ia][1]*e3[1] + e3_refs[ts][ia][2]*e3[2]);

					//} catch(char const* s) {
					} else {
						std::cout << "ERROR: ts,ia " << ts << " " << ia << " " << e1_flag << " " << e2_flag << " " << e3_flag << std::endl;
						e1_dtd0.add(1.0);
						e2_dtd0.add(1.0);
						e3_dtd0.add(1.0);
					}
				}
			}

			const double e1x  = e1_dtd0.value;
			const double e2x  = e2_dtd0.value;
			const double e3x  = e3_dtd0.value;
	
			e1_corr[local_time].add(e1x, e1_dtd0.num_samples);
			e2_corr[local_time].add(e2x, e2_dtd0.num_samples);
			e3_corr[local_time].add(e3x, e3_dtd0.num_samples);

			e1_l2_corr[local_time].add(0.5*(3.0*e1x*e1x-1.0), e1_dtd0.num_samples);
			e2_l2_corr[local_time].add(0.5*(3.0*e2x*e2x-1.0), e2_dtd0.num_samples);
			e3_l2_corr[local_time].add(0.5*(3.0*e3x*e3x-1.0), e3_dtd0.num_samples);

			/*
			std::cout << "num_samples,e1,e2,e3: " << 
				e1_corr[local_time].num_samples << " " <<  
				e2_corr[local_time].num_samples << " " <<  
				e3_corr[local_time].num_samples << std::endl;
			*/
			
		}

	}

	void save(std::string filename="rot_corr.dat")
	{
		//std::cout << "saving corr data. e1_corr.size() " <<  e1_corr.size() << std::endl;
		
		std::ofstream fout_corr(filename);	
		fout_corr << "time   e1_corr   e2_corr   e3_corr   e1_l2_corr   e2_l2_corr   e3_l2_corr   e1_samples   e2_samples   e3_samples" << std::endl;
		std::cout << "time   e1_corr   e2_corr   e3_corr   e1_l2_corr   e2_l2_corr   e3_l2_corr   e1_samples   e2_samples   e3_samples" << std::endl;
		for(int ts = 0; ts <  e1_corr.size(); ts++)
		{
			double e1 = 0.0, e2 = 0.0, e3 = 0.0; 
			double e1_l2 = 0.0, e2_l2 = 0.0, e3_l2 = 0.0; 

			if (e1_corr[ts].num_samples > 0) e1 = e1_corr[ts].value/e1_corr[ts].num_samples; 
			if (e2_corr[ts].num_samples > 0) e2 = e2_corr[ts].value/e2_corr[ts].num_samples; 
			if (e3_corr[ts].num_samples > 0) e3 = e3_corr[ts].value/e3_corr[ts].num_samples; 

			if (e1_corr[ts].num_samples > 0) e1_l2 = e1_l2_corr[ts].value/e1_l2_corr[0].value;
			if (e2_corr[ts].num_samples > 0) e2_l2 = e2_l2_corr[ts].value/e2_l2_corr[0].value; 
			if (e3_corr[ts].num_samples > 0) e3_l2 = e3_l2_corr[ts].value/e3_l2_corr[0].value; 

			fout_corr << ts*dt_per_frame << " " << e1 << " " << e2 << " " << e3 << " " <<
				e1_l2 << " " << e2_l2 << " " << e3_l2 << " " <<
			       	e1_corr[ts].num_samples << " " << 
			       	e2_corr[ts].num_samples << " " << 
			       	e3_corr[ts].num_samples << std::endl;
			std::cout << ts*dt_per_frame << " " << e1 << " " << e2 << " " << e3 << " " <<
				e1_l2 << " " << e2_l2 << " " << e3_l2 << " " <<
			       	e1_corr[ts].num_samples << " " << 
			       	e2_corr[ts].num_samples << " " << 
			       	e3_corr[ts].num_samples << std::endl;

		}

		fout_corr.close();
	}

};

auto angle_N1HN2 = make_histogram(axis::regular<double>(180, 0.0, 180.0));
auto angle_N1HN3 = make_histogram(axis::regular<double>(180, 0.0, 180.0));

class hbond_stats
{
	std::string datadir;
	int corr_length, corr_interval;
	float dt_per_frame; 

	int total_tseries, total_measurement_frames; 

	std::vector<fs::path> filepaths; 

	NeighborList first_nbr; 

	bond_stat bstat; 

	std::string file_histogram, file_correlation, file_hbond_pairs;

        // num_alive_bonds, ratio, d0dt, p1, p2
	typedef std::tuple<int, double, double, double, double> corr_t;
	std::map<int,corr_t> corr;

	hbond_populations hbpop; 
	hb_profiler hb_prof;

	rotational_correlations rot_corr; 

public:

	hbond_stats(
			std::string _datadir, 
			int _corr_length, 
			int _corr_interval, 
			float _dt_per_frame, 
			std::string file_corr, 
			std::string file_histo, 
			std::string file_hpair) : 

		datadir(_datadir), 
		corr_length(_corr_length), 
		corr_interval(_corr_interval), 
		dt_per_frame(_dt_per_frame), 
		file_correlation(file_corr), 
		file_histogram(file_histo), 
		file_hbond_pairs(file_hpair)
	{
		auto cpath = fs::current_path();
		std::cout << "cpath : " << cpath.string() << "\n";

		std::vector<fs::path> filepaths_all;

		// read filelist from file or directory
		if(fs::is_regular_file(datadir))
		{
			std::ifstream fin(datadir);
			std::string line;

			while(std::getline(fin,line))
				filepaths_all.push_back(line);

		} else {
			for (auto f : fs::recursive_directory_iterator(datadir))
			{
				if(f.path().string().find(".xyz") != std::string::npos)
				{
					filepaths_all.push_back(f.path());
					//std::cout << f.path().string() << "\n"; 
				}
			}
		}


		std::sort(begin(filepaths_all),end(filepaths_all), 
			[](auto const &p1, auto const &p2) {return p1.string() < p2.string();} );

		for(int i = 0; i < filepaths_all.size(); i++)
			filepaths.push_back(filepaths_all[i]);

		// get initial frame
		std::string filename = filepaths[0].string();
		std::ifstream ifile(filename);	
		MDFrame single_frame = read_single_mdframe(ifile, filename);
		first_nbr = NeighborList(single_frame);

		total_measurement_frames = corr_length + ((filepaths.size() - corr_length)/corr_interval)*corr_interval;
		total_tseries = 1 + (filepaths.size() - corr_length)/corr_interval;

		std::cout << filename <<  "\n=============================================================\n";
		std::cout << "total_frames, total_tseries, total_measurement_frames : " 
			<< filepaths.size() << " " << total_tseries << " " << total_measurement_frames << std::endl;
		std::cout << filename <<  "=============================================================\n";

		//hbpop.initialize(filepaths.size(), single_frame.natoms);
		hbpop.initialize(corr_length, single_frame.natoms);

//		hb_prof.initialize(single_frame.natoms, filepaths);

		rot_corr.initialize(corr_length, corr_interval, dt_per_frame, total_tseries);
	}

	void calc_corr()
	{
		std::ofstream fout_correlation(file_correlation);	

		for (int iframe = 0; iframe <filepaths.size(); iframe++)
		{
			//if(iframe >= total_measurement_frames) continue; 
			if(iframe > total_measurement_frames) continue; 
			if(iframe%corr_interval==0) hbpop.allocate_new_hbond_map_all(); 
	
			std::string filename = filepaths[iframe].string();
			std::ifstream ifile(filename);	
			MDFrame single_frame = read_single_mdframe(ifile, filename);
			//single_frame.print();
	
			auto nbr = NeighborList(single_frame);
			//nbr.print();
	
			assert(first_nbr.nbrlist_ncenter.size() == nbr.nbrlist_ncenter.size()); 

			std::cout << filename <<  " ===========================\n";

			// register new reference time 
			if(iframe%corr_interval == 0) rot_corr.register_ref_vectors(nbr.nbrlist_ncenter);

			rot_corr.compute_correlations(iframe, nbr.nbrlist_ncenter);

		
//			hb_prof.get_neighbor_n_distances(iframe, nbr);
//			hb_prof.get_average_hb_distance(iframe, nbr);

			for (int tseries_id = 0; tseries_id < hbpop.get_num_tseries(); tseries_id++)
			{
				int local_step  = iframe - tseries_id*corr_interval; 
				assert(local_step >= 0);

				if(local_step >= corr_length) continue; 

				/*
				std::cout << "tseries_id,local_step,iframe,clength,cinterval: " 
					<< tseries_id << " " <<  local_step << " " <<  iframe << " "
					<< corr_length << " " << corr_interval << std::endl;
				*/

				//=== distribution calculations ===//
				for (int iatom = 0; iatom< nbr.nbrlist_hcenter.size(); iatom++)
				{
					///std::cout << "i,nbrlist_hcenter.size() " << iatom <<  " " << nbr.nbrlist_hcenter[iatom].size() << std::endl;
	
					const auto n = nbr.nbrlist_hcenter[iatom]; 
	
					/*
					// N1-H-N3 angle
					if(n.size() >= 3) 
					{
						int j = 0;  // N1
						int jgid = std::get<1>(n[j]);
						const auto nj_r = std::get<3>(n[j]);
						const double nj_norm = std::sqrt(nj_r.dx*nj_r.dx + nj_r.dy*nj_r.dy + nj_r.dz*nj_r.dz);
	
						int k = 2;  // N3
						int kgid = std::get<1>(n[k]);
						const auto nk_r = std::get<3>(n[k]);
						const double nk_norm = std::sqrt(nk_r.dx*nk_r.dx + nk_r.dy*nk_r.dy + nk_r.dz*nk_r.dz);
		
						double dot = nj_r.dx*nk_r.dx + nj_r.dy*nk_r.dy + nj_r.dz*nk_r.dz;
						double cosine = dot/nj_norm/nk_norm;
						double degree = std::acos(cosine) * 180.0 / M_PI;
	
						if(tseries_id==0) angle_N1HN3(degree); // once in the tseries loop.
	
						//std::cout << "jgid,kgid,nj_norm,nk_norm,N3_s.C2: " << jgid << " " << kgid << " " << nj_norm << " " << nk_norm << " " << 
						//	hbpop.N3_s.C2 << " " << (nk_norm < hbpop.N3_s.C2) << std::endl;
	
						if(nk_norm < hbpop.N3_s.C2) {
							hbpop.N3_s.sample(tseries_id, local_step, jgid, kgid, nk_norm);
							//if(hbpop.N2_s2l.C1 < nk_norm) hbpop.N3_s.sample(tseries_id, local_step, jgid, kgid);
						}
					}
					*/
	
					// N1-H-N2 angle
					if(n.size() >= 2) 
					{
						int j = 0;  // N1
						int jgid = std::get<1>(n[j]);
						const auto nj_r = std::get<3>(n[j]);
						const double nj_norm = std::sqrt(nj_r.dx*nj_r.dx + nj_r.dy*nj_r.dy + nj_r.dz*nj_r.dz);
	
						int k = 1;  // N2 
						int kgid = std::get<1>(n[k]);
						const auto nk_r = std::get<3>(n[k]);
						const double nk_norm = std::sqrt(nk_r.dx*nk_r.dx + nk_r.dy*nk_r.dy + nk_r.dz*nk_r.dz);
		
						double dot = nj_r.dx*nk_r.dx + nj_r.dy*nk_r.dy + nj_r.dz*nk_r.dz;
						double cosine = dot/nj_norm/nk_norm;
						double degree = std::acos(cosine) * 180.0 / M_PI;
	
						if(tseries_id==0) angle_N1HN2(degree); // once in the tseries loop
	
						if(nk_norm < hbpop.N2_s.C2) hbpop.N2_s.sample(tseries_id, local_step, jgid, kgid, nk_norm);

						//std::cout << "jgid,kgid,nj_norm,nk_norm,N2_s.C2: " << jgid << " " << kgid << " " << nj_norm << " " << nk_norm << " " << 
						//	hbpop.N2_s.C2 << " " << (nk_norm < hbpop.N2_s.C2) << std::endl;
					}
				}

			}
	
		}
		rot_corr.save();

		fout_correlation.close();
//		hb_prof.save();

	};

	void save_result()
	{
		std::vector<std::vector<double>> angles;

		// normalize all histograms
		auto all_histograms = {angle_N1HN2, angle_N1HN3};

		std::cout << "num samples "; 
		for(auto h : all_histograms)
			std::cout << std::accumulate(h.begin(), h.end(), 0.0);
		std::cout << std::endl;

		for(auto h : all_histograms)
		{
			std::vector<double> angle;

			double sum = std::accumulate(h.begin(), h.end(), 0.0);
			for(int i=0; i<h.axis().size(); i++) angle.push_back(h[i]/sum);
			angles.push_back(angle);
		}

		std::ofstream fout_histogram(file_histogram);	
	
		fout_histogram << "angle    N1-H-N2    N1-H-N3 " << std::endl;
		for(int idx=0; idx<angles[0].size(); idx++) 
		{
			fout_histogram << idx << " ";
			for(int j=0; j<angles.size(); j++) fout_histogram << angles[j][idx] << " ";
			fout_histogram << std::endl;
		}

		fout_histogram.close();

		std::ofstream fout_correlation(file_correlation);	

		std::vector<std::vector<double>> corrs;
		std::cout << "N2_s:\n"; corrs.push_back(hbpop.N2_s.summary());
		//std::cout << "N3_s:\n"; corrs.push_back(hbpop.N3_s.summary());

		std::ofstream fout_hbond_pairs(file_hbond_pairs);
		//fout_hbond_pairs << "Time(fs)     Prob(N2)      Prop(N3)     Samples(N2)     Samples(N3) " << std::endl;
		fout_hbond_pairs << "Time(fs)     Prob(N2)      Samples(N2)     " << std::endl;

		for(int i = 0; i < corrs[0].size(); i++) 
			fout_hbond_pairs << i*dt_per_frame << " " << 
			corrs[0][i]/corrs[0][0] << " " << corrs[0][i] << std::endl;
			//corrs[0][i]/corrs[0][0] << " " << corrs[1][i]/corrs[1][0] << " " << corrs[0][i] << " " << corrs[1][i] << " " << std::endl;

		fout_hbond_pairs.close();

	};

};

int main(int argc, char* argv[])
{

	if(argc != 5)
	{
		std::cout << "USAGE: ./a.out datadir corr_length corr_interval dt_per_frame" << std::endl;
		exit(0);
	}

	std::string datadir(argv[1]);
	float corr_length = std::atof(argv[2]);
	float corr_interval = std::atof(argv[3]);
	float dt_per_frame = std::atof(argv[4]);

	int corr_length_frame = int(corr_length/dt_per_frame);
	int corr_interval_frame = int(corr_interval/dt_per_frame);

	std::cout << "=============================================================" << std::endl;
	std::cout << "datadir : " << datadir << std::endl; 
	std::cout << "corr_length(fs), corr_length(frame)  : " << corr_length << " " << corr_length_frame << std::endl; 
	std::cout << "corr_interval(fs), corr_interval(frame) : " << corr_interval << " " << corr_interval_frame << std::endl; 
	std::cout << "dt_per_frame(fs/frame) : " << dt_per_frame << std::endl; 
	std::cout << "=============================================================" << std::endl;

	hbond_stats hstats(datadir, corr_length_frame, corr_interval_frame, dt_per_frame, "corr.dat", "hbond.dat", "hbpop.dat");
	hstats.calc_corr();
	hstats.save_result();

	return 0;
}
