
// Keep in private GIT
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <thread>
#include <deque>
#include <unordered_set>
using namespace std;

#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FParameterNames.hpp"
enum TaskType {P2M = 0, M2M, M2L, M2L_OUT, L2L, L2P, P2P, P2P_OUT};
string taskNames[] = {"P2M", "M2M", "M2L", "M2L_out", "L2L", "L2P", "P2P", "P2P_out"};
struct Task
{
	TaskType type;
	long long int uniqueId;
	vector<long long int> id;
	int mpiNode;
	int level;
	bool operator==(const Task & other) const
	{
		if(type != other.type || id.size() != other.id.size())
			return false;
		for(size_t i = 0; i < id.size(); ++i)
			if(id[i] != other.id[i])
				return false;
		return true;
	}
	bool operator!=(const Task & other) const
	{
		return !((*this)==other);
	}
	void print(void)
	{
		cout << taskNames[type];
		for(size_t i = 0; i < id.size(); ++i)
			cout << ", " << id[i];
		cout << endl;
	}
};
//Spécialisation de hash pour le type Task
namespace std {
  template <> struct hash<Task>
  {
    size_t operator()(const Task & x) const
    {
		return x.uniqueId;
    }
  };
}
struct DagData
{
	unordered_set<Task> allTask;
};

bool parseLine(DagData & dagData, deque<string> & lineElements)
{
	if(lineElements.size() < 1)
		return false;
	Task task;
	if(lineElements.size() >= 14 && lineElements[0] == "P2P_out")
	{
		task.type = P2P_OUT;
		task.uniqueId = stoll(lineElements[1]);
		task.id.resize(4);
		task.id[0] = stoll(lineElements[9]);
		task.id[1] = stoll(lineElements[10]);
		task.id[2] = stoll(lineElements[11]);
		task.id[3] = stoll(lineElements[12]);
		task.mpiNode = stoi(lineElements[13]);
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 10 && lineElements[0] == "P2P")
	{
		task.type = P2P;
		task.uniqueId = stoll(lineElements[1]);
		task.id.resize(4);
		task.id[0] = stoll(lineElements[5]);
		task.id[1] = stoll(lineElements[6]);
		task.id[2] = stoll(lineElements[7]);
		task.id[3] = stoll(lineElements[8]);
		if(task.id[0] == 0 && task.id[1] == 0 && task.id[2] == 0 && task.id[3] == 0)
			cout << "Suricate" << endl;
		task.mpiNode = stoi(lineElements[9]);
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 10 && lineElements[0] == "M2L" )
	{
		task.type = M2L;
		task.uniqueId = stoll(lineElements[1]);
		task.id.resize(5);
		task.id[0] = stoll(lineElements[2]);
		task.id[1] = stoll(lineElements[5]);
		task.id[2] = stoll(lineElements[6]);
		task.id[3] = stoll(lineElements[7]);
		task.id[4] = stoll(lineElements[8]);
		task.mpiNode = stoi(lineElements[9]);
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 13 && lineElements[0] == "M2L_out")
	{
		task.type = M2L_OUT;
		task.uniqueId = stoll(lineElements[1]);
		task.id.resize(5);
		task.id[0] = stoll(lineElements[2]);
		task.id[1] = stoll(lineElements[8]);
		task.id[2] = stoll(lineElements[9]);
		task.id[3] = stoll(lineElements[10]);
		task.id[4] = stoll(lineElements[11]);
		task.mpiNode = stoi(lineElements[12]);
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 13 && lineElements[0] == "M2M")
	{
		task.type = M2M;
		task.uniqueId = stoll(lineElements[1]);
		task.id.resize(5);
		task.id[0] = stoll(lineElements[2]);
		task.id[1] = stoll(lineElements[8]);
		task.id[2] = stoll(lineElements[9]);
		task.id[3] = stoll(lineElements[10]);
		task.id[4] = stoll(lineElements[11]);
		task.mpiNode = stoi(lineElements[12]);
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 13 && lineElements[0] == "L2L")
	{
		task.type = L2L;
		task.uniqueId = stoll(lineElements[1]);
		task.id.resize(5);
		task.id[0] = stoll(lineElements[2]);
		task.id[1] = stoll(lineElements[8]);
		task.id[2] = stoll(lineElements[9]);
		task.id[3] = stoll(lineElements[10]);
		task.id[4] = stoll(lineElements[11]);
		task.mpiNode = stoi(lineElements[12]);
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 8 && lineElements[0] == "L2P")
	{
		task.type = L2P;
		task.uniqueId = stoll(lineElements[1]);
		task.id.resize(2);
		task.id[0] = stoll(lineElements[5]);
		task.id[1] = stoll(lineElements[6]);
		task.mpiNode = stoi(lineElements[7]);
		dagData.allTask.insert(task);
	}
	else if(lineElements.size() >= 8 && lineElements[0] == "P2M")
	{
		task.type = L2P;
		task.uniqueId = stoll(lineElements[1]);
		task.id.resize(2);
		task.id[0] = stoll(lineElements[5]);
		task.id[1] = stoll(lineElements[6]);
		task.mpiNode = stoi(lineElements[7]);
		dagData.allTask.insert(task);
	}
	else
	{
		cout << "No match for " << lineElements[0] << " - " << lineElements.size() << endl;
		return false;
	}
	return true;
}
void split(string& line, string const& delim, deque<string> &result)
{
	size_t prevPos = 0;
	size_t pos = 0;
	while((pos = line.find(delim, prevPos)) != string::npos)
	{
		result.push_back(line.substr(prevPos, pos-prevPos));
		prevPos = pos+delim.size();
	}
	if(prevPos < line.size())
		result.push_back(line.substr(prevPos));
}
bool fillDagData(const char* const filename, DagData & dagData)
{
	std::ifstream fichier(filename, ios::in);  // on ouvre le fichier en lecture

	if(!fichier)  // si l'ouverture a réussi
	{       
		cerr << "Couldn't open " << filename << endl;
		return false;
	}
	string line;
	string delim(", ");
	deque<string> splitLine;
	int count = 0;
	while(!fichier.eof())
	{
		++count;
		getline(fichier, line);
		splitLine.clear();
		split(line, delim, splitLine);
		parseLine(dagData, splitLine);
	}
	cout << (count-1) << " lines in " << filename << endl;
	// instructions
	fichier.close();  // on ferme le fichier
	return true;
}
void compareDag(DagData& dag1, DagData& dag2, int treeHeight)
{
	long long int notFoundCount[treeHeight] = {0};
	bool notFound[treeHeight] = {false};
	for(Task task : dag1.allTask)
	{
		bool found = false;
		if(task.type == P2P || task.type == P2P_OUT || task.type == P2M || task.type == L2P)
			notFound[treeHeight-1] = true;
		else if(task.id[0] < treeHeight)
			++notFound[task.id[0]] = true;
		for(Task task2 : dag2.allTask)
		{
			if(task == task2)
			{
				found = true;
				break;
			}
		}
		if(found == false)
		{
			task.print();
			if(task.type == P2P || task.type == P2P_OUT || task.type == P2M || task.type == L2P)
				++notFoundCount[treeHeight-1];
			else
				++notFoundCount[task.id[0]];
		}
	}
	for(int i = 0; i < treeHeight; ++i)
		if(notFound[i] == true)
			cout << "Diff lvl " << i << " -> " << notFoundCount[i] << endl;
}
int main(int argc, char* argv[]){
    const FParameterNames Explicit {
        {"-e"},
        "Trace from explicit mpi"
    };
	const FParameterNames Implicit {
		{"-i"} ,
		"Trace from implicit mpi"
	};
    const FParameterNames TreeHeight {
        {"-h"},
        "Height of the tree"
    };
    FHelpDescribeAndExit(argc, argv, "Compare DAG mapping", Explicit, Implicit, TreeHeight);

    // Get params
    const char* const explicitFilename = FParameters::getStr(argc,argv,Explicit.options, "explicit.rec");
    const char* const implicitFilename = FParameters::getStr(argc,argv,Implicit.options, "implicit.rec");
    const int treeHeight = FParameters::getValue(argc,argv,TreeHeight.options, 5);

	DagData implicitData, explicitData;
	bool implicitGood, explicitGood;
	std::thread explicitThread([&](){
		explicitGood = fillDagData(explicitFilename, explicitData);
		});
	explicitThread.join();
	std::thread implicitThread([&](){
		implicitGood = fillDagData(implicitFilename, implicitData);
		});
	implicitThread.join();
	if(implicitGood && explicitGood)
	{
		cout << explicitData.allTask.size() << " tasks in explicit." << endl;
		cout << implicitData.allTask.size() << " tasks in implicit." << endl;
		compareDag(explicitData, implicitData, treeHeight);
	}
    return 0;
}
