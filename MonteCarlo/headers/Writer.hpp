#ifndef ____Writer____
#define ____Writer____

// File Stream and etc.
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

// Error catch
#include <signal.h>
#include <cstdlib>

#ifdef _WIN32
#include <windows.h>
#endif // _WIN32
#ifdef linux
#include <sys/stat.h>
#endif // linux




using namespace std;

class Writer {
    public:
        ofstream myfile;
        string file_name;
        string new_file_name;
        string directory_name;
        void MakeDirectory();
        // Wirter(string filename);
        // ~Writer();
        Writer(string filename, int filenum, string directory_name);
        // Writer(bool Foo){/*Intentialnally empty for FooWriter*/};
        void FindNextFileNum(int filenum);
        void OpenNewFile();
        template<typename T>
        void WriteLine(vector<T> contents);
        void WriteLine(vector<string> contents);
        void WriteLine(string contents);
        void WriteLine(double contents);
        void WriteLine(long double contents);
        void WriteLine(int contents);
        void CloseNewFile();
        static vector<double> Argument_reader(string file_path,int size){
            ifstream mf(file_path);
            if(mf.fail()){
                cerr << "Not correct file path" << endl;
                exit(100);
            }
            vector<double> res;
            for(int i = 0; i < size; i++){
                string var_name, equal_sign;
                mf >> var_name >> equal_sign; // Not use.
                double data;
                mf >> data;
                res.push_back(data);
            }
            mf.close();
            return res;
        };
        static vector<long double> Argument_readerL(string file_path,int size){
            ifstream mf(file_path);
            if(mf.fail()){
                cerr << "Not correct file path" << endl;
                exit(100);
            }
            vector<long double> res;
            for(int i = 0; i < size; i++){
                string var_name, equal_sign;
                mf >> var_name >> equal_sign; // Not use.
                long double data;
                mf >> data;
                res.push_back(data);
            }
            mf.close();
            return res;
        };
};

// Writer::Writer(string filename){
//     this-> new_file_name = filename;
//     this-> directory_name = "Result";
// }

Writer::Writer(string filename, int filenum = 1, string directory_name = "Result"){
    this-> file_name = filename;
    this-> directory_name = directory_name;

    this-> MakeDirectory();
    this-> FindNextFileNum(filenum);
    this-> OpenNewFile();
    cout << "Saving Start: " << this-> new_file_name << "\n";
    myfile.precision(10);
}
// Writer::~Writer(){
//     this-> CloseNewFile();
// }

void Writer::MakeDirectory(){
    #ifdef _WIN32
    CreateDirectory(this->directory_name.c_str(), NULL);
    #endif // _WIN32

    #ifdef linux
    mkdir(this->directory_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    #endif // linux
}

void Writer::FindNextFileNum(int filenum = 1){
    ifstream f(this-> new_file_name);
    int cf = filenum;
    while(true){
        this-> new_file_name = this->file_name + "_" + to_string(cf++) + ".csv";
        ifstream f(this-> new_file_name);
        if(!f.good()) break;
    }
    f.close();
}

void Writer::OpenNewFile(){
    this->myfile.open(this->new_file_name);
}
template<typename T>
void Writer::WriteLine(vector<T> contents){
    // This assumes that there is no line break.
    myfile << fixed << setw(3) << int(contents.at(0)) << ',';
    for(int i = 1; i < contents.size()-1; i++)
        myfile << fixed << setw(9) << contents.at(i) << ',';
    myfile << fixed << setw(9) << contents.back() << '\n';
}

void Writer::WriteLine(vector<string> contents){
    // This assumes that there is no line break.
    myfile << fixed << setw(3) << contents.at(0) << ',';
    for(int i = 1; i < contents.size()-1; i++)
        myfile << fixed << setw(11) << contents.at(i) << ',';
    myfile << fixed << setw(11) << contents.back() << '\n';
}

void Writer::WriteLine(string contents){
    // This assumes that there is no line break.
    myfile << contents;
}

void Writer::WriteLine(long double contents){
    // This assumes that there is no line break.
    myfile << contents;
}
void Writer::WriteLine(double contents){
    // This assumes that there is no line break.
    myfile << contents;
}

void Writer::WriteLine(int contents){
    // This assumes that there is no line break.
    myfile << contents;
}

void Writer::CloseNewFile(){
    this-> myfile.close();
    cout << "Save Completed: " << this-> new_file_name << "\n";
}

class FooWriter : Writer {
    public:
        ofstream myfile;
        string file_name;
        string new_file_name;
        string directory_name;
        void MakeDirectory(){}
        // FooWriter(string filename, string directory_name = "Result") : Writer(true) {}
        FooWriter(string filename, int filenum = 1, string directory_name = "Result") : Writer(filename,filenum,directory_name) {}
        void FindNextFileNum(){}
        void OpenNewFile(){}
        void WriteLine(string contents){}
        void WriteLine(double contents){}
        void WriteLine(int contents){}
        void CloseNewFile(){}
};

#endif // ____Writer____
