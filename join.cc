
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
#include <string.h>
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
std::map<int,std::map<int,double> > read_graph(char* file)
{
	FILE* fp;
	char line[4092];
	std::string summary;
	std::string tmp;
	int v_max;
	int e_max;
	int a;
	int b;
	double c;

	fp = fopen (file , "r");
	summary = fgets(line, sizeof(line), fp);
	std::stringstream ss;

	ss << summary;
	ss >> tmp >> tmp >> v_max >> e_max;
	
	std::map<int,std::map<int,double> > G;
	std::stringstream entry;

	for (int i=0; i<v_max; i++)
	{
		G[i] = std::map<int,double>({});
	}
	int i = 0;
	while (i<e_max)
	{
		fgets(line, sizeof(line), fp);
		if ((char)line[0] == 'c')
		{
			continue;
		}
		ss << line;
		ss >> tmp >> a >> b >> c;
		G[a-1].insert(std::pair<int,double>(b-1,c));
		G[b-1].insert(std::pair<int,double>(a-1,c));
		i++;
	}
	fclose(fp);
	return G;
}

//read in the Graph (non - dimacs) (set the base index to 0 instead of 1)
std::map<int,std::map<int,double> > read_graph_no_dimacs(char const * file, std::map<int,std::map<int,double>> Gold, std::map<int,int> IC)
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
	
	std::map<int,std::map<int,double> > G;
	std::stringstream entry;

	for (int i=0; i<v_max; i++)
	{
		G[i] = std::map<int,double>({});
	}

	for (int i=0; i<e_max; i++)
	{
		fgets(line, sizeof(line), fp);
		ss << line;
		ss >> a >> b;
		G[IC.at(a+1)].insert(std::pair<int,double>(IC.at(b+1),Gold.at(IC.at(a+1)).at(IC.at(b+1))));
		G[IC.at(b+1)].insert(std::pair<int,double>(IC.at(a+1),Gold.at(IC.at(b+1)).at(IC.at(a+1))));
	}
	fclose(fp);
	return G;
}


//print graph in dimacs format (fix base index to 1 again)
void dmx2tmp(std::map<int,std::map<int,double>> G, std::string tmp, std::map<int,int> &IC)
{
	std::string s1 = "";
	int count = 0;
	int iter = 0;
	std::map<int,int> out {};
	for (auto const v : G)
	{
		if (IC.find(v.first) == IC.end())
		{
			iter++;
			IC[v.first]=iter;
			out[iter]=v.first;
		}
		for (auto const w : v.second)
		{
			if (v.first > w.first)
			{
				continue;
			}
			if (IC.find(w.first) == IC.end())
			{
				iter++;
				IC[w.first]=iter;
				out[iter]=w.first;
			}
			s1.append("\ne ");
			s1.append(to_string(IC[v.first]));
			s1.append(" ");
			s1.append(to_string(IC[w.first]));
			s1.append(" ");
			s1.append(to_string(int(w.second*10000.0)));
			//s1.append(to_string(w.second));
			count++;
		}
	}
	std::string s2 = "p edge ";
	s2.append(to_string(IC.size()));
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
	IC = out;
}

//print graph in dimacs format (fix base index to 1 again)
void dmx2tmp2(std::map<int,std::map<int,double>> G, std::string tmp)
{
	std::string s1 = "";
	int count = 0;
	for (auto const v : G)
	{
		for (auto const w : v.second)
		{
			if (v.first < w.first)
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

//find edge in stack
int FindInStack(std::map<int,std::map<int,double>> stack, std::pair<int,int> edge)
{
	if (stack.find(edge.first)==stack.end())	
	{
		if (stack[edge.first].find(edge.second)==stack[edge.first].end())
		{
			return 0;
		}
	}
	return 1;
}

//symmetric difference on edge sets
std::map<int,std::map<int,double>> E_symmetric(std::vector<std::map<int,std::map<int,double>>> vG)
{
	std::map<int,std::map<int,double>> good = vG.at(0);
	std::map<int,std::map<int,double>> bad  = {};
	dmx2tmp2(good,"./good1.dmx");
	for ( unsigned int i=1; i<vG.size(); i++)
	{
		std::cout<<"HEEEEEEEEELLLLLLLLLLLOOOOOOOOOOOOO"<<std::endl;
		for (auto v : vG.at(i))
		{
			for ( auto w : v.second )
			{
				if (FindInStack(bad,std::pair<int,int>(v.first,w.first))==1)
				{
					continue;
				}
				if (FindInStack(good,std::pair<int,int>(v.first,w.first))==1)
				{
					std::cout<<"out 1"<<std::endl;
					good[v.first].erase(w.first);
				}
				else if (good.find(v.first) != good.end())
				{
					std::cout<<"in 1"<<std::endl;
					good[v.first].insert(w);
				}
				else
				{
					std::cout<<"in 2"<<std::endl;
					good[v.first] = {};
					good[v.first].insert(w);
				}
			}
		}
	}
	for ( auto v : good )
	{
		if (v.second.size() == 0)
		{
			good.erase(v.first);
		}
	}
	dmx2tmp2(good,"./good.dmx");
	return good;
}

//reconstruct shortest path
std::map<int,std::map<int,double>> reconstruct(int a, int b, std::map<int,std::map<int,int>> P, std::map<int,std::map<int,double>> G)
{
	std::cout<<"bam"<<std::endl;
	std::map<int,std::map<int,double>> out = {};
	int i=b;
	while (true)
	{
		std::cout<<"bam2"<<std::endl;
		if (P[a][i] == a)
		{
			out[i].insert(std::pair<int,double>(a,G[i][a]));
			out[a].insert(std::pair<int,double>(i,G[a][i]));
			std::cout<<"bam3"<<std::endl;
			dmx2tmp2(out,"./duh.dmx");
			return out;
		}
		out[i].insert(std::pair<int,double>(a,G[i][P[a][i]]));
		out[P[a][i]].insert(std::pair<int,double>(P[a][i],G[P[a][i]][i]));
		i = P[a][i];
	}
}

//find symetric difference of shortest paths of perfect matching
std::map<int,std::map<int,double>> SymPath(std::map<int,std::map<int,double>> PM, std::map<int,std::map<int,int>> P, std::map<int,std::map<int,double>> G)
{
	std::map<int,std::map<int,double>> out = {};
	std::vector<std::map<int,std::map<int,double>>> tmp = {}; 
	for (auto const v : PM)
	{
		for (auto const w : v.second)
		{
			std::cout<<"bam1"<<std::endl;
			tmp.push_back(reconstruct(v.first,w.first,P,G));
		}
	}
	return E_symmetric(tmp);
}

//find weigh and length of t-join
void findWeight(std::map<int,std::map<int,double>> G, double &weight, int &length)
{
	for (auto v : G)
	{
		for (auto w : v.second)
		{
			if (v.first < w.first)
			{
				weight+=w.second;
				length++;
			}
		}
	}
}

//find vertices with negative weight edges
void Gpos_Tminus(std::map<int,std::map<int,double>> G, std::map<int,std::map<int,double>> &Gpos, std::vector<int> &Tm,std::map<int,std::map<int,double>> &Em)
{
	int count;
	std::map<int,double> tmp;
	std::map<int,double> tmp2;
	for (auto const v : G)
	{
		tmp = {};
		tmp2 = {};
		count =1;
		for (auto const w : v.second)
		{
			if (0.0>w.second)
			{
				count*=-1;
				tmp.insert(std::pair<int,double>(w.first,-w.second));
				if (v.first < w.first)
				{
					tmp2.insert(std::pair<int,double>(w.first,w.second));
				}
			}
			else
			{
				tmp.insert(w);
			}
		}
		Gpos.insert(std::pair<int,std::map<int,double>>(v.first,tmp));
		Em.insert(std::pair<int,std::map<int,double>>(v.first,tmp2));
		std::cout<<count<<std::endl;
		if (count == -1)
		{
			Tm.push_back(v.first);
		}
	}
}

//reduction step in algorithm
std::map<int,std::map<int,double>> reduce(std::map<int,std::map<int,double>> G)
{
	std::map<int,std::map<int,double>> out;
	std::map<int,double> tmp;
	double m=G.begin()->second.begin()->second;
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
			tmp.insert(std::pair<int,double>(w.first,w.second - m));
		}
		out.insert(std::pair<int,std::map<int,double>>(v.first,tmp));
	}
	return out;
}

// add weight of t-join to all edges
void reduce_by(std::map<int,std::map<int,double>> &G, double C)
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
std::pair<int,std::map<int,double>> dijk(std::map<int,std::map<int,double> > G, int s, std::map<int,int> &P, std::vector<int> T)
{
	std::map<int,double> l;
	std::vector<int> R = {};
	std::vector<int> V = {};
	int min;

	for (unsigned int i=0;i<G.size();i++)
	{
		V.push_back(i);
		l.insert(std::pair<int,double>(i,std::numeric_limits<double>::max()));
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
					P[iter.at(i)] = min;
				}
			}
		}
	}
	l.erase(s);
	std::map<int,double> L;
	for (auto const v : l)
	{
		if (std::find(T.begin(), T.end(), v.first) == T.end())
		{
			continue;
		}
		if (std::numeric_limits<int>::max() != v.second)
		{
			L.insert(v);
		}
	}
	return std::pair<int,std::map<int,double>> {s,L};
}

//metric closure
std::map<int,std::map<int,double>> metric_closure(std::map<int,std::map<int,double> > G, std::vector<int> T, std::map<int,std::map<int,int>> &P)
{
	std::map<int,std::map<int,double>> Gbar = {};
	std::pair<int,std::map<int,double>> D;
	std::map<int,int> tmp;
	for (auto v : G)
	{
		if (std::find(T.begin(),T.end(), v.first) == T.end())
		{
			continue;
		}
		tmp = {};
		std::cout<<"hello b"<<std::endl;
		Gbar.insert(dijk(G,v.first,tmp,T));
		P[v.first] = tmp;
		std::cout<<"hello c"<<std::endl;
		P.insert(std::pair<int,std::map<int,int>>(v.first,tmp));
	}
	return Gbar;		
}

//minimum weight perfect matching
std::map<int,std::map<int,double>> MWPMP(std::map<int,std::map<int,double>> G)
{
	std::map<int,int> IC {};
	//write_file
	dmx2tmp2(G,"./surprise.dmx");
	dmx2tmp(G,"./tmp.dmx",IC);
	std::cout<<"hello 1.5"<<std::endl;
	//do system call to blossom5
	system("rm ./tmp2.dmx 2> /dev/null; touch ./tmp2.dmx");
	system("./blossom5 -V -e ./tmp.dmx -w ./tmp2.dmx >./junk.txt");
	std::cout<<"hello2"<<std::endl;
	//read file back in
	std::map<int,std::map<int,double>> M = read_graph_no_dimacs("./tmp2.dmx",G,IC);
	dmx2tmp2(M,"./dump.dmx");
	std::cout<<"hello2"<<std::endl;
	//system("rm ./tmp* 2> /dev/null");
	return M;
}

//find minimum weight empty set join
std::map<int,std::map<int,double>> min_tjoin(std::map<int,std::map<int,double> > G)
{
	//initialize parameters
	std::vector<int> Tminus = {};
	std::map<int,std::map<int,double>> Eminus = {};
	std::map<int,std::map<int,double>> Gpos = {};
	std::map<int,std::map<int,int>> P = {};
	std::map<int,std::map<int,double>> out = {};
	std::vector<std::map<int,std::map<int,double>>> tmp = {};

	//apply values
	Gpos_Tminus(G,Gpos,Tminus,Eminus);
	dmx2tmp2(Eminus,"./eminus.dmx");
	dmx2tmp2(Gpos,"./gpos.dmx");

	std::cout<<"hello a"<<std::endl;
	std::map<int,std::map<int,double>> Gbar = metric_closure(Gpos,Tminus,P);
	std::cout<<"hello b"<<std::endl;
	print_vector(Tminus);
	std::map<int,std::map<int,double>> MWPM = MWPMP(Gbar);
	dmx2tmp2(MWPM,"./pow.dmx");
	//symmetric difference on edge sets
	tmp.push_back(SymPath(MWPM,P,G));
	dmx2tmp2(SymPath(MWPM,P,Gpos),"./trying.dmx");
	tmp.push_back(Eminus); 
	std::cout<<"now!"<<std::endl;
	out = E_symmetric(tmp); 
	return out;
}

//main algorithm
void main_algorithm(std::map<int,std::map<int,double> > G)
{
	std::map<int,std::map<int,double>> ejoin;
	std::map<int,std::map<int,double>> rG = reduce(G);
	dmx2tmp2(rG,"./blah2.dmx");
	double weight;
	int length;
	int counter=0;
	while (counter<1)
	//while (true)
	{
		//compute min t-join
		weight = 0.0;
		length = 0;
		std::cout<<"hello1"<<std::endl;
		dmx2tmp2(rG,"./blah.dmx");
		ejoin = min_tjoin(rG);
		findWeight(ejoin,weight,length);	
		dmx2tmp2(ejoin,"./joinold.dmx");
		std::cout<<"hello n"<<std::endl;
		counter++;
		
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
				reduce_by(rG, -weight/((double)length));
			}
		}
		else
		{
			reduce_by(rG, -weight/((double)length));
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
		std::map<int,std::map<int,double> > G = read_graph(argv[1]);
		//run algorithm
		main_algorithm(G);
		return 0;
	}
}

// The End //

