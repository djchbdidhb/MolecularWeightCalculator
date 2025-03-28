#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

using std::cin;
using std::cout;
using std::endl;
using std::setw;
using std::left;
using std::right;
using std::map;
using std::vector;
using std::string;
using std::pair;
using std::istream;
using std::find_if;
using std::find;
using std::sort;
using std::isdigit;
using std::istringstream;
using std::ifstream;
using std::invalid_argument;
using std::runtime_error;

inline bool find_element_begin(const char c)
{
    return (c >= 'A' && c<= 'Z') || ( c == '(' || c == '[' || c == '{');
}
inline bool if_not_digit(const char c)
{
    return !isdigit(c);
}
inline bool if_not_lowercase(const char c)
{
    return !(c >= 'a' && c<= 'z');
}

map<string, int> read_input(istream&);
double calculate_mass(const map<string, int>&, const map<string, double>&);
void convert_symbol(map<string, int>&);
string::const_iterator find_bracket_end(string::const_iterator, string::const_iterator);
map<string, int> deal_bracket(const string&, string::const_iterator);
void print_result(const map<string, int>&, const map<string, double>&, double);
map<string, double> read_atom_mass(void);
int string_to_int(const string&);

int main()
{
    map<string, double> aw;
    try{
        aw = read_atom_mass();
    }catch(runtime_error e){
        cout << e.what() << endl;
        return 0;
    }
    cout << "Enter the chemical formula: ";
    map<string, int> atom_num = read_input(cin);
    try{
        convert_symbol(atom_num);
        double mass = calculate_mass(atom_num, aw);
        print_result(atom_num, aw, mass);
    }catch(invalid_argument e1){
        cout << e1.what() << endl;
    }catch(runtime_error e2){
        cout << e2.what() << endl;
    }
    return 0;
}

void convert_symbol(map<string, int>& atom_num)
{
    ifstream input("symbol.txt");
    if(!input){
        throw runtime_error("Cannot open symbol.txt");
    }
    string line;
    string abbrev;
    while(getline(input, line)){
        istringstream iss(line);
        iss >> abbrev;
        if(atom_num.find(abbrev) != atom_num.end()){
            string formula;
            iss >> formula;
            istringstream f_iss(formula);
            map<string, int> temp = read_input(f_iss);
            for(map<string, int>::iterator j = temp.begin();
            j != temp.end(); j++){
                j->second *= atom_num[abbrev];
                atom_num[j->first] += j->second;
            }
            atom_num.erase(abbrev);
        }
    }
}

void print_result(const map<string, int>& atom_num,
    const map<string, double>& aw, double mass)
{
    vector<pair<string, int> > elements(atom_num.begin(), atom_num.end());
    sort(elements.begin(), elements.end(),
    [&aw](const pair<string, int>& a, const pair<string, int>& b){
        return aw.at(a.first) < aw.at(b.first);
    });
    cout << left << setw(12) << "Element";
    cout << left << setw(12) << "Count";
    cout << "Atomic Weight" << endl;
    cout << string(40, '-') << endl;
    for(vector<pair<string, int> >::const_iterator i = elements.begin();
    i != elements.end(); i++){
        cout << left << setw(12) << i->first;
        cout << left << setw(12) << i->second;
        cout << aw.at(i->first) << endl;
    }
    cout << string(40, '-') << endl;
    cout << "Molecular Mass: " << mass << endl;
}

double calculate_mass(const map<string, int>& m, const map<string, double>& aw)
{
    double total_mass = 0;
    for(map<string, int>::const_iterator i = m.begin(); i != m.end(); i++){
        if(aw.find(i->first) == aw.end()){
            throw invalid_argument("Unsupported symbol: " + i->first);
        }
        total_mass += aw.find(i->first)->second * i->second;
    }
    return total_mass;
}

map<string, int> read_input(istream& in)
{
    map<string, int> ret;
    string s;
    in >> s;
    string::const_iterator i = s.begin();
    string::const_iterator j;
    while(i != s.end()){
        if(*i == '(' || *i == '[' || *i == '{'){
            j = find_bracket_end(i + 1, s.cend());
            istringstream iss(string(i + 1, j));
            map<string, int> temp = read_input(iss);
            j++;
            if(j != s.end() && isdigit(*j)){
                string::const_iterator k = find_if(j, s.cend(), if_not_digit);
                int x = string_to_int(string(j, k));
                for(map<string, int>::iterator it = temp.begin();
                it != temp.end(); it++){
                    it->second *= x;
                }
                j = k;
            }
            for(map<string, int>::const_iterator it = temp.begin();
            it != temp.end(); it++){
                ret[it->first] += it->second;
            }
            i = j;
        }else if(i != s.end()){
            j = find_if(i+1, s.cend(), if_not_lowercase);
            if(isdigit(*j)){
                string::const_iterator k = find_if(j, s.cend(), if_not_digit);
                ret[string(i, j)] += string_to_int(string(j, k));
                i = j = k;
            }else{
                ++ret[string(i, j)];
                i = j;
            }
        }
        i = find_if(i, s.cend(), find_element_begin);
    }
    return ret;
}

map<string, double> read_atom_mass(void)
{
    map<string, double> aw;
    ifstream input("atomic_weight.txt");
    if(!input){
        throw runtime_error("Cannot open atomic_weight.txt");
    }
    string line;
    while(getline(input, line)){
        istringstream iss(line);
        string element_name;
        string temp;
        iss >> temp >> element_name >> temp;
        iss >> aw[element_name];
    }
    return aw;
}

string::const_iterator find_bracket_end(string::const_iterator begin, 
    string::const_iterator end)
{
    int num_of_bracket = 1;
    string::const_iterator k = begin;
    for(; k != end; k++){
        if(*k == '(' || *k == '[' || *k == '{'){
            num_of_bracket++;
        }else if(*k == ')' || *k == ']' || *k == '}'){
            num_of_bracket--;
        }
        if(num_of_bracket == 0){
            break;
        }
    }
    return k;
}

int string_to_int(const string& s)
{
    istringstream iss(s);
    int x;
    iss >> x;
    return x;
}