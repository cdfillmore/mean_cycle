
//////////////////////////////////////////////////////
//                                                  //
//  Minimum Mean Weight Cycle Algorithm             //
//                                                  //
//  written by Christopher Fillmore                 //
//             Giridhar Shenoy                      //
//             Nadezhda Vassilyeva                  //
//                                                  //
//  to compile: g++ join.cc -o join -std=c++11      //
//  to execute: join w_test1.dmx		    //
//                                                  //
//////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <map>
#include <algorithm>

using namespace std;

//intersection of vectors
std::vector<int> intersection(std::vector<int> v1, std::vector<int> v2)
{
	std::vector<int> v3;
	
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	
	std::set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(v3));
	
	return v3;
}

//union of vectors
std::vector<int> union_1(std::vector<int> v1, std::vector<int> v2)
{
	std::vector<int> v3;
	
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	
	std::set_union(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(v3));
	
	return v3;
}

//symmetric difference of vectors
std::vector<int> symmetric(std::vector<int> v1, std::vector<int> v2)
{
	std::vector<int> v3;
	
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	
	std::set_symmetric_difference(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(v3));
	
	return v3;
}

//print a vector
void print_vector(std::vector<int> a)
{
	std::string out1="";
	for (unsigned int i=0; i<a.size(); i++)
	{
		out1.append(" ");
		out1.append(std::to_string(a.at(i)));
	}
	std::cout<<out1<<std::endl;
}

//read in the Graph (dimacs) (set the base index to 0 instead of 1)
std::map<int,std::map<int,float> > read_graph(char* file)
{
	FILE* fp;
	char line[4092];
	std::string summary;
	std::string tmp;
	int v_max;
	int e_max;
	int a;
	int b;
	float c;

	fp = fopen (file , "r");
	summary = fgets(line, sizeof(line), fp);
	std::stringstream ss;

	ss << summary;
	ss >> tmp >> tmp >> v_max >> e_max;
	
	std::map<int,std::map<int,float> > G;
	std::stringstream entry;

	for (int i=0; i<v_max; i++)
	{
		G[i] = std::map<int,float>({});
	}

	for (int i=0; i<e_max; i++)
	{
		fgets(line, sizeof(line), fp);
		ss << line;
		ss >> tmp >> a >> b >> c;
		G[a-1].insert(std::pair<int,float>(b-1,c));
		G[b-1].insert(std::pair<int,float>(a-1,c));
	}
	fclose(fp);
	return G;
}

//read in the Graph (non - dimacs) (set the base index to 0 instead of 1)
std::map<int,std::map<int,float> > read_graph_no_dimacs(char const * file, std::map<int,std::map<int,float>> Gold)
{
	FILE* fp;
	char line[4092];
	std::string summary;
	std::string tmp;
	int v_max;
	int e_max;
	int a;
	int b;

	fp = fopen (file , "r");
	summary = fgets(line, sizeof(line), fp);
	std::stringstream ss;

	ss << summary;
	ss >> v_max >> e_max;
	
	std::map<int,std::map<int,float> > G;
	std::stringstream entry;

	for (int i=0; i<v_max; i++)
	{
		G[i] = std::map<int,float>({});
	}

	for (int i=0; i<e_max; i++)
	{
		fgets(line, sizeof(line), fp);
		ss << line;
		ss >> a >> b;
		G[a-1].insert(std::pair<int,float>(b,Gold.at(a).at(b)));
		G[b-1].insert(std::pair<int,float>(a,Gold.at(b).at(a)));
	}
	fclose(fp);
	return G;
}


//print graph in dimacs format (fix base index to 1 again)
void dmx2tmp(std::map<int,std::map<int,float>> G, std::vector<int> T, std::string tmp)
{
	std::string s1 = "";
	int count = 0;
	for (auto const v : G)
	{
		if (std::find(T.begin(),T.end(),v.first) != T.end())
		{
			for (auto const w : v.second)
			{
				if (std::find(T.begin(),T.end(),w.first) != T.end())
				{
					s1.append("\ne ");
					s1.append(to_string(v.first+1));
					s1.append(" ");
					s1.append(to_string(w.first+1));
					s1.append(" ");
					s1.append(to_string(w.second));
					count++;
				}
			}
		}
	}
	std::string s2 = "p edge ";
	s2.append(to_string(T.size()));
	s2.append(" ");
	s2.append(to_string(count));
	s2.append(s1);
	
	std::string cmd = "rm ";
	cmd.append(tmp);
	cmd.append(" 2> /dev/null; touch ");
	cmd.append(tmp);
	system(cmd.c_str());

	ofstream myfile;
	myfile.open(tmp.c_str());
	myfile << s2;
	myfile.close();
}

//print graph in dimacs format (fix base index to 1 again)
void dmx2tmp2(std::map<int,std::map<int,float>> G, std::string tmp)
{
	std::string s1 = "";
	int count = 0;
	for (auto const v : G)
	{
		for (auto const w : v.second)
		{
			s1.append("\ne ");
			s1.append(to_string(v.first+1));
			s1.append(" ");
			s1.append(to_string(w.first+1));
			s1.append(" ");
			s1.append(to_string(w.second));
			count++;
		}
	}
	std::string s2 = "p edge ";
	s2.append(to_string(G.size()));
	s2.append(" ");
	s2.append(to_string(count));
	s2.append(s1);
	
	std::string cmd = "rm ";
	cmd.append(tmp);
	cmd.append(" 2> /dev/null; touch ");
	cmd.append(tmp);
	system(cmd.c_str());

	ofstream myfile;
	myfile.open(tmp.c_str());
	myfile << s2;
	myfile.close();
}


//symmetric difference on edge sets
std::map<int,std::map<int,float>> E_symmetric(std::map<int,std::map<int,float>> G1, std::map<int,std::map<int,float>> G2, float &weight, int &length)
{
	std::map<int,std::map<int,float>> out = {};
	std::map<int,float> tmp1;
	std::map<int,float> tmp2;
	float SUM;
	for (auto const v : G1)
	{
		tmp1 = {};
		if (G2.find(v.first) == G2.end())
		{
			tmp1 = v.second;
			SUM = 0.0;
			for(std::map<int,float>::iterator it = tmp1.begin(); it != tmp1.end(); ++it) SUM += it->second;
			weight += SUM;
			length += tmp1.size();
		}
		else
		{
			for (auto const w : v.second)
			{
				if (G2.at(v.first).find(w.second) == G2.at(v.first).end())
				{
					tmp1.insert(w);
					weight += w.second;
					length++;
				}
			}
		}
		out.insert(std::pair<int,std::map<int,float>>(v.first,tmp1));
	}
	for (auto const v : G2)
	{
		if (G1.find(v.first) == G1.end())
		{
			tmp2 = v.second;
			SUM = 0.0;
			for(std::map<int,float>::iterator it = tmp2.begin(); it != tmp2.end(); ++it) SUM += it->second;
			weight += SUM;
			length += tmp2.size();
			
		}
		else
		{
			for (auto const w : v.second)
			{
				if (G1.at(v.first).find(w.second) == G1.at(v.first).end())
				{
					tmp2.insert(w);
					weight += w.second;
					length++;
				}
			}
		}
		out.insert(std::pair<int,std::map<int,float>>(v.first,tmp2));
	}
	return out;
}

//find vertices with negative weight edges
void Gpos_Tminus(std::map<int,std::map<int,float>> G, std::map<int,std::map<int,float>> &Gpos, std::vector<int> &Tm,std::map<int,std::map<int,float>> &Em)
{
	int count;
	std::map<int,float> tmp;
	std::map<int,float> tmp2;
	for (auto const v : G)
	{
		tmp = {};
		count =0;
		for (auto const w : v.second)
		{
			if (0.0>w.second)
			{
				count = 1;
				tmp.insert(std::pair<int,float>(w.first,-w.second));
				tmp2.insert(std::pair<int,float>(w.first,w.second));
			}
			else
			{
				tmp.insert(w);
			}
		}
		Gpos.insert(std::pair<int,std::map<int,float>>(v.first,tmp));
		Em.insert(std::pair<int,std::map<int,float>>(v.first,tmp2));
		if (count == 1)
		{
			Tm.push_back(v.first);
		}
	}
}

//reduction step in algorithm
std::map<int,std::map<int,float>> reduce(std::map<int,std::map<int,float>> G)
{
	std::map<int,std::map<int,float>> out;
	std::map<int,float> tmp;
	float m=G.begin()->second.begin()->second;
	for (auto const v : G)
	{
		for (auto const w : v.second)
		{
			if (m<w.second)
			{
				m = w.second;
			}
		}
	}
	for (auto v : G)
	{
		tmp = {};
		for (auto w : v.second)
		{
			tmp.insert(std::pair<int,float>(w.first,w.second - m));
		}
		out.insert(std::pair<int,std::map<int,float>>(v.first,tmp));
	}
	return out;
}

// add weight of t-join to all edges
void reduce_by(std::map<int,std::map<int,float>> &G, float C)
{
	for (auto v : G)
	{
		for (auto w : v.second)
		{
			G.at(v.first).at(w.first) += C;
		}
	}
}

//dijkstra's algorithm
std::pair<int,std::map<int,float>> dijk(std::map<int,std::map<int,float> > G, int s)
{
	std::map<int,float> l;
	std::vector<int> R = {};
	std::vector<int> V = {};
	int min;

	for (unsigned int i=0;i<G.size();i++)
	{
		V.push_back(i);
		l.insert(std::pair<int,float>(i,std::numeric_limits<float>::max()));
	}
	l[s] = 0.0;
	std::vector<int> iter = intersection(V,symmetric(R,V));
	while (iter.size() > 0 )
	{
		min = iter.at(0);
		for (unsigned int i=1;i<iter.size();i++)
		{
			if (l.at(iter.at(i)) < l.at(min))
			{
				min=iter.at(i);
			}
		}
		R = union_1(R,std::vector<int> {min});
		iter = intersection(V,symmetric(R,V));
		for (unsigned int i=0;i<iter.size();i++)
		{
			if (G.at(min).find(iter.at(i))!=G.at(min).end())
			{
				if (l.at(iter.at(i)) > (l.at(min) + G.at(min).at(iter.at(i))))
				{
					l.at(iter.at(i)) = l.at(min) + G.at(min).at(iter.at(i));
				}
			}
		}
	}
	l.erase(s);
	std::map<int,float> L;
	for (auto const v : l)
	{
		if (std::numeric_limits<int>::max() != v.second)
		{
			L.insert(v);
		}
	}
	return std::pair<int,std::map<int,float>> {s,L};
}

//metric closure
std::map<int,std::map<int,float>> metric_closure(std::map<int,std::map<int,float> > G)
{
	std::map<int,std::map<int,float>> Gbar = {};
	for (unsigned int i=0;i<G.size();i++)
	{
		Gbar.insert(dijk(G,i));
	}
	return Gbar;		
}

//minimum weight perfect matching
std::map<int,std::map<int,float>> MWPMP(std::map<int,std::map<int,float>> G, std::vector<int> T)
{
	print_vector(T);
	//write_file
	dmx2tmp(G,T,"./tmp.dmx");
	//do system call to blossom5
	system("rm ./tmp2.dmx 2> /dev/null; touch ./tmp2.dmx");
	system("./blossom5 -V -e ./tmp.dmx -w ./tmp2.dmx >./junk.txt");
	std::cout<<"hello2"<<std::endl;
	//read file back in
	std::map<int,std::map<int,float>> M = read_graph_no_dimacs("./tmp2.dmx",G);
	std::cout<<"hello2"<<std::endl;
	system("rm ./tmp* 2> /dev/null");
	return M;
}

//find minimum weight empty set join
std::map<int,std::map<int,float>> min_tjoin(std::map<int,std::map<int,float> > G, float &weight, int &length)
{
	//initialize parameters
	std::vector<int> T = {};
	std::vector<int> Tminus = {};
	std::map<int,std::map<int,float>> Eminus = {};
	std::map<int,std::map<int,float>> Gpos = {};

	//apply values
	Gpos_Tminus(G,Gpos,Tminus,Eminus);
	std::vector<int> T2 = symmetric(T,Tminus);

	std::map<int,std::map<int,float>> Gbar = metric_closure(Gpos);
	std::map<int,std::map<int,float>> MWPM = MWPMP(Gbar,T2);
	//symmetric difference on edge sets
	return E_symmetric(MWPM,Eminus,weight,length); 
}

//main algorithm
void main_algorithm(std::map<int,std::map<int,float> > G)
{
	std::map<int,std::map<int,float>> ejoin;
	std::map<int,std::map<int,float>> rG = reduce(G);
	dmx2tmp2(rG,"./blah2.dmx");
	float weight;
	int length;
	int counter=0;
	while (true)
	{
		//compute min t-join
		weight = 0.0;
		length = 0;
		std::cout<<"hello1"<<std::endl;
		ejoin = min_tjoin(rG,weight,length);
		std::cout<<"hello n"<<std::endl;
		counter++;
		
		dmx2tmp2(rG,"./blah.dmx");
		std::cout<<weight<<std::endl;
		if ( weight < 1e-5 )
		{
			if (weight > -1e-5)
			{
				dmx2tmp2(ejoin,"./join.dmx");
				break;
			}
			else
			{
				dmx2tmp2(ejoin,"./joinold.dmx");
				reduce_by(rG, -weight/((float)length));
			}
		}
		else
		{
			reduce_by(rG, -weight/((float)length));
		}
	}
	
}

//main function
int main(int argc, char* argv[])
{
	//check args for dmx file
	if (argc < 1)
	{
		cout << "no graph specified" << endl;
	return -1;
	}
	else
	{
		system("echo \"hi\"");
		//read Graph
		std::map<int,std::map<int,float> > G = read_graph(argv[1]);
		//run algorithm
		main_algorithm(G);
		return 0;
	}
}

// The End //

