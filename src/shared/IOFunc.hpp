#ifndef JDSURFG_SHARED_IOFUNC_H_
#define JDSURFG_SHARED_IOFUNC_H_

#include <fstream>
#include <regex.h>
#include "numerical.hpp"

void create_directory(const char *dirname);
void print_progressbar(float percentage);


/**
 * @brief template function to read parameters by using regex
 * @example read_par_regex("NX",nx,"Par_file")
 * 
 * @tparam T datatype
 * @param filename filename
 * @param varname  variable name
 * @param var      reference to a variable
 */
template<typename T> int
read_par_regex(const std::string &varname,T &var,std::ifstream &infile)
{   
    // go to beginning of the file
    infile.clear();
    infile.seekg(0);

    // temporary vars
    std::string tempname = "^" + varname,dummy;
    std::istringstream info;
    int ierr = 1;

    // read line by line and find the match str
    std::string line;
    while(getline(infile,line)){
        if(line.length() == 0 || line[0] == '#') continue;;

        // match 
        regex_t reg;
        regmatch_t pmatch[1];
        regcomp(&reg,tempname.data(),REG_NEWLINE);
        int status = regexec(&reg,line.c_str(),1,pmatch,0);
        if(status == 0){
            info.str(line);
            info >> dummy; info >> dummy;
            info >> var;
            info.clear();
            ierr = 0;
            break;
        }
        regfree(&reg);
    }

    if(ierr == 1){
      printf("cannot find %s\n",varname.c_str());
    }
    
    return ierr;
};



/**
 * @brief template function to read string parameters by using regex
 * @example read_par_regex("NX",nx,"Par_file")
 * 
 * @tparam T datatype
 * @param filename filename
 * @param varname  variable name
 * @param var      reference to a variable
 */
template<> inline int 
read_par_regex<std::string>(const std::string &varname,std::string &var,
                            std::ifstream &infile)
{   
    // go to beginning of the file
    infile.clear();
    infile.seekg(0);

    // temporary vars
    std::string tempname = "^" + varname;
    int ierr = 1;

    // read line by line and find the match str
    std::string line;
    while(getline(infile,line)){
        if(line.length() == 0 || line[0] == '#') continue;;

        // match 
        regex_t reg;
        regmatch_t pmatch[1];
        regcomp(&reg,tempname.data(),REG_NEWLINE);
        int status = regexec(&reg,line.c_str(),1,pmatch,0);
        if(status == 0){
            size_t i = 0;
            for(; i < line.size(); i ++ ){
                if(line[i] == '=') break;
            }
            var = line.substr(i+1);
            ierr = 0;
            break;
        }
        regfree(&reg);
    }

    if(ierr == 1){
      printf("cannot find %s\n",varname.c_str());
    }
    
    return ierr;
}

#endif // end JDSURFG_SHARED_IOFUNC_H_
