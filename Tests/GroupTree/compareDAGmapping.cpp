
// Keep in private GIT
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <thread>
using namespace std;

#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FParameterNames.hpp"
struct DagData
{

};
void parseLine(DagData & dagData, string& line)
{
	std::regex regexAct("^([a-z0-9]+)_"); 
	std::regex regexArgP2M("^[a-z0-9]+_([0-9]+)_([0-9]+)_([0-9])$"); 
	std::regex regexArgM2M("^[a-z0-9]+_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)$"); 
	std::regex regexArgP2P("^[a-z0-9]+_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)$"); 
	std::smatch matches;
	if(std::regex_search (line,matches,regexAct))
	{
		if(matches[1] == "p2m")
		{
			std::regex_search (line,matches,regexArgP2M);
			int startMorton = stoi(matches[1], nullptr,10);
			int endMorton = stoi(matches[2], nullptr,10);
			int node = stoi(matches[3], nullptr,10);
		}
		else if(matches[1] == "m2m")
		{
			std::regex_search (line,matches,regexArgM2M);
			int level = stoi(matches[1], nullptr,10);
			int startMorton = stoi(matches[2], nullptr,10);
			int endMorton = stoi(matches[3], nullptr,10);
			int startMorton2 = stoi(matches[4], nullptr,10);
			int endMorton2 = stoi(matches[5], nullptr,10);
			int node = stoi(matches[6], nullptr,10);
		}
		else if(matches[1] == "m2l")
		{
			std::regex_search (line,matches,regexArgM2M);
			int level = stoi(matches[1], nullptr,10);
			int startMorton = stoi(matches[2], nullptr,10);
			int endMorton = stoi(matches[3], nullptr,10);
			int startMorton2 = stoi(matches[4], nullptr,10);
			int endMorton2 = stoi(matches[5], nullptr,10);
			int node = stoi(matches[6], nullptr,10);
		}
		else if(matches[1] == "l2l")
		{
			std::regex_search (line,matches,regexArgM2M);
			int level = stoi(matches[1], nullptr,10);
			int startMorton = stoi(matches[2], nullptr,10);
			int endMorton = stoi(matches[3], nullptr,10);
			int startMorton2 = stoi(matches[4], nullptr,10);
			int endMorton2 = stoi(matches[5], nullptr,10);
			int node = stoi(matches[6], nullptr,10);
		}
		else if(matches[1] == "l2p")
		{
			std::regex_search (line,matches,regexArgP2M);
			int startMorton = stoi(matches[1], nullptr,10);
			int endMorton = stoi(matches[2], nullptr,10);
			int node = stoi(matches[3], nullptr,10);
		}
		else if(matches[1] == "p2p")
		{
			std::regex_search (line,matches,regexArgP2P);
			int startMorton = stoi(matches[1], nullptr,10);
			int endMorton = stoi(matches[2], nullptr,10);
			int startMorton2 = stoi(matches[3], nullptr,10);
			int endMorton2 = stoi(matches[4], nullptr,10);
			int node = stoi(matches[5], nullptr,10);
		}
	}
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
	while(!fichier.eof())
	{
		getline(fichier, line);
		if(line.size() > 3 && line[0] == 'N')
		{
			line = line.substr(3);
			parseLine(dagData, line);
		}
	}
	// instructions
	fichier.close();  // on ferme le fichier
	return true;
}
void compareDag(DagData& dag1, DagData& dag2)
{
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
    FHelpDescribeAndExit(argc, argv, "Compare DAG mapping", Explicit, Implicit);

    // Get params
    const char* const explicitFilename = FParameters::getStr(argc,argv,Explicit.options, "explicit.rec");
    const char* const implicitFilename = FParameters::getStr(argc,argv,Implicit.options, "implicit.rec");

	DagData implicitData, explicitData;
	bool implicitGood, explicitGood;
	std::thread implicitThread([&](){
		explicitGood = fillDagData(explicitFilename, explicitData);
		});
	std::thread explicitThread([&](){
		implicitGood = fillDagData(implicitFilename, implicitData);
		});
	implicitThread.join();
	explicitThread.join();
	if(implicitGood && explicitGood)
		compareDag(implicitData, explicitData);
    return 0;
}
