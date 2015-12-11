
//////////////////////////////////////////////////////
//                                                  //
//  Edmonds max. cardinality matching algorithm     //
//                                                  //
//  written by Christopher Fillmore                 //
//             Giridhar Shenoy                      //
//             Nadezhda Vassilyeva                  //
//                                                  //
//  to compile: g++ Blossom.cc -o bloss -std=c++11  //
//  to execute: bloss test.dmx                      //
//                                                  //
//////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
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

//read in the Graph (set the base index to 0 instead of 1)
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

//print the matching in dimacs format (fix base index to 1 again)
void print_vector(std::vector<int> V)
{
	std::string tmp;
	int v_max = V.size();
	int e_max = 0;
	std::vector<int> visited = {};
	for (int i=0; i<v_max; i++)
	{
		if (std::find(visited.begin(),visited.end(),i) == visited.end())
		{
			e_max++;
			visited.push_back(V.at(i));
		}
	}
	std::cout<<"p edge "<<v_max<<e_max<<std::endl;
	for (int i=0; i<e_max; i++)
	{
		if (visited.at(i) > V.at(visited.at(i)))
		{
			std::cout<<"e "<<V.at(visited.at(i))+1<<" "<<visited.at(i)+1<<std::endl;
		}
		else
		{
			std::cout<<"e "<<visited.at(i)+1<<" "<<V.at(visited.at(i))+1<<std::endl;
		}
	}
	std::cout<<tmp<<std::endl;
}

void reduce(std::map<int,std::map<int,float>>& G)
{
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
		for (auto w : v.second)
		{
			G.at(v.first).at(w.first) = w.second - m;
		}
	}
}

//main algorithm
void main_algorithm(std::map<int,std::map<int,float> > G, int mean)
{
	//initialize parameters
	int h=0;
	reduce(G);
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
		//set mean to zero
		int mean = 0;
		main_algorithm(G,mean);
		//run algorithm
		//std::vector<int> mu = find_max_mu(G,mean);
		//print the mu vector
		//std::cout<<mean<<std::endl;
		//print_vector(mu);
		return 0;
	}
}

// The End //

