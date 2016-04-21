
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <cassert>
#include <sstream>
#include <cstring>

template <class VariableType>
inline VariableType StrToOther(const std::string& str, const VariableType& defaultValue = VariableType(), bool* hasWorked = nullptr){
    std::istringstream iss(str,std::istringstream::in);
    VariableType value = defaultValue;
    iss >> value;
    if(hasWorked){
        if(!iss.eof()){
            char c;
            iss >> c;
            (*hasWorked) = (c == '\n' || c == '\0');
        }
        else{
            (*hasWorked) = true;
        }
    }
    if( /*iss.tellg()*/ iss.eof() ) return value;
    return defaultValue;
}


std::string ReduceName(const std::string name){
    const std::vector<std::pair<std::string,std::string>> mapping=
    { {"M2M-level", "M2M"} , {"M2L-level", "M2L"} , {"M2L-out-level", "M2L-out"} ,
      {"L2L-level", "L2L"} };

    for(const std::pair<std::string,std::string> mp : mapping){
        if(name.substr(0, mp.first.length()) == mp.first){
            return mp.second;
        }
    }
    return name;
}

struct LineData{
    std::string name;
    int nb;
    std::string type;
    double duration;

    LineData(const char line[]){
        std::vector<std::string> words;
        {
            int start = 0;
            int end = 1;
            while(line[end] != '\0'){
                while(line[end] != '\0'
                      && line[end] != ','){
                    end += 1;
                }

                if(line[end] != '\0'){
                    words.push_back(std::string(&line[start], end-start));
                    end += 1;
                    start = end;
                }
            }
            if(start != end){
                words.push_back(std::string(&line[start], end-start));
            }
        }
        if(words.size() != 4){
            printf("Error line is no composed of 4 words\n");
            exit(111);
        }
        name = ReduceName(words[0].substr(1, words[0].size() - 2));
        bool hasWorked = false;
        nb = StrToOther<int>(words[1], -1, &hasWorked);
        if(hasWorked == false){
            printf("Error cannot get nb val\n");
            exit(112);
        }
        type = words[2].substr(1, words[2].size() - 2);        
        duration = StrToOther<double>(words[3], 0, &hasWorked);
        if(hasWorked == false){
            printf("Error cannot get duration val\n");
            exit(112);
        }
    }
};

struct TimeData{
    TimeData()
        : tt(0), tr(0), ti(0){
    }

    double tt; // time in task
    double tr; // time in runtime
    double ti; // time idle
};

/*
"Name","Count","Type","Duration"
"Initializing",24,"Runtime",1.259594
"Overhead",6304,"Runtime",150.667381
"Idle",1635,"Other",2.964436
"Scheduling",3285,"Runtime",50.173307
"Sleeping",330,"Other",10978.895357
"FetchingInput",1611,"Runtime",2.834788
"execute_on_all_wrapper",48,"Task",60.98058
"PushingOutput",1611,"Runtime",49.857017
"P2P-out",403,"Task",3840.070231
"Callback",1447,"Runtime",0.459004
"P2M",49,"Task",7759.555381
"M2M-level-5",49,"Task",916.867961
"M2L-level-5",7,"Task",3866.40312
"M2M-level-4",7,"Task",90.183757
"M2L-out-level-5",32,"Task",809.783766
"P2P",49,"Task",28378.015095
"L2L-level-5",49,"Task",749.115965
"M2L-level-6",49,"Task",33069.582498
"M2L-out-level-6",806,"Task",10014.65321
"L2P",49,"Task",10532.198512
"Deinitializing",24,"Runtime",0.600177
"M2L-level-3",1,"Task",45.115451
"M2L-level-2",1,"Task",2.638928
"L2L-level-2",1,"Task",1.343462
"L2L-level-4",7,"Task",87.756298
"L2L-level-3",1,"Task",10.658414
"M2M-level-3",1,"Task",11.480571
"M2M-level-2",1,"Task",1.41104
"M2L-level-4",1,"Task",511.345345
  */

int main(int argc, char** argv){
    if(argc != 4){
        printf("Error usage is:\n"
               "%s seq_file parallel_file%%d nb_threads\n",
               argv[0]);
        return 200;
    }

    printf("seq file is %s\n", argv[1]);
    printf("parallel file are %s\n", argv[2]);
    const int nbThreads = StrToOther<int>(argv[3], -1);
    if(nbThreads == 1){
        printf("Error cannot convert nb threads\n");
        return 201;
    }
    printf("up to %d threads\n", nbThreads);


    std::vector<std::unordered_map<std::string,double>> timeTasks;
    std::unordered_set<std::string> allTaskNames;
    std::vector<TimeData> times;

    for(int idxFile = 0 ; idxFile <= nbThreads ; ++idxFile){
        char filename[1024];
        if(idxFile == 0){
            strncpy(filename, argv[1], 1024);
        }
        else{
            sprintf(filename, argv[2], idxFile);
        }
        timeTasks.emplace_back();
        times.emplace_back();

        printf("Open file : %s\n", filename);
        FILE* timeFile = fopen(filename, "r");
        if(timeFile == NULL){
            printf("Cannot open file %s\n", filename);
            return 99;
        }
        size_t sizeLine = 1024;
        char* line = (char*)malloc(sizeLine);

        {// Read header
            if((sizeLine = getline((char**)&line, &sizeLine, timeFile)) == -1){
                printf("Cannot read header\n");
                return 1;
            }
            // Should be: "Name","Count","Type","Duration"
            if(strcmp("\"Name\",\"Count\",\"Type\",\"Duration\"\n",
                      line) != 0){
                printf("Header is incorrect\n");
                return 2;
            }
        }

        while((sizeLine = getline((char**)&line, &sizeLine, timeFile)) != -1){
            LineData dt(line);
            // Task, Runtime, Other
            if(dt.type == "Task"){
                if(dt.name != "execute_on_all_wrapper"){
                    timeTasks[idxFile][dt.name] += dt.duration;
                    allTaskNames.insert(dt.name);
                    times[idxFile].tt += dt.duration;
                }
            }
            else if(dt.type == "Runtime"){
                if(dt.name == "Scheduling"
                        || dt.name == "FetchingInput"
                        || dt.name == "PushingOutput"){
                    times[idxFile].tr += dt.duration;
                }
            }
            else if(dt.type == "Other"){
                if(dt.name == "Idle"){
                    times[idxFile].ti += dt.duration;
                }
            }
            else {
                printf("Arg do not know type %s\n", dt.type.c_str());
                return 3;
            }
        }

        fclose(timeFile);
    }

    // Global efficiencies
    {
        // Manually set seq idel and runtime
        times[0].ti = 0;
        times[0].tr = 0;
        times[1].ti = 0;
        times[1].tr = 0;

        printf("Create global-eff.data\n");

        FILE* resFile = fopen("global-eff.data", "w");
        assert(resFile);
        fprintf(resFile, "0 \tgranularity-eff \ttasks-eff \truntime-eff \tpipeline-eff\n");

        for(int idx = 1; idx <= nbThreads ; ++idx){
            fprintf(resFile, "%d \t%e \t%e \t%e \t%e\n",
                    idx,
                    times[0].tt/times[idx].tt,
                    times[1].tt/times[idx].tt,
                    times[idx].tt/(times[idx].tt+times[idx].tr),
                    (times[idx].tt+times[idx].tr)/(times[idx].tt+times[idx].tr+times[idx].ti));
        }

        fclose(resFile);
    }

    // Global efficiencies
    {
        printf("Create task-eff.data\n");

        FILE* resFile = fopen("task-eff.data", "w");
        assert(resFile);
        fprintf(resFile, "0 ");
        for(const std::string tsk : allTaskNames){
            fprintf(resFile, "\t%s ", tsk.c_str());
        }
        fprintf(resFile, "\n");

        for(int idx = 1; idx <= nbThreads ; ++idx){
            fprintf(resFile, "%d", idx);


            for(const std::string tsk : allTaskNames){
                fprintf(resFile, "\t%e ", timeTasks[1][tsk]/timeTasks[idx][tsk]);
            }
            fprintf(resFile, "\n");
        }

        fclose(resFile);
    }

    {
        printf("Create task-gr-eff.data\n");

        FILE* resFile = fopen("task-gr-eff.data", "w");
        assert(resFile);
        fprintf(resFile, "0 ");
        for(const std::string tsk : allTaskNames){
            fprintf(resFile, "\t%s ", tsk.c_str());
        }
        fprintf(resFile, "\n");

        for(int idx = 1; idx <= nbThreads ; ++idx){
            fprintf(resFile, "%d", idx);


            for(const std::string tsk : allTaskNames){
                fprintf(resFile, "\t%e ", timeTasks[0][tsk]/timeTasks[idx][tsk]);
            }
            fprintf(resFile, "\n");
        }

        fclose(resFile);
    }

    return 0;
}
