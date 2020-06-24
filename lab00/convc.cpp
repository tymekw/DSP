/*
// okolo 1000 przy uzyciu vector

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

int main(){

    vector<int> x(24), h(24);
    ifstream f_read("x_h_sign.txt");

    for(int i=0; i <24; ++i){
        f_read >> x[i] >> h[i];
    }

    /*
    for(int i=0; i < x.size(); i++)
      std::cout << x.at(i) << ' ';


    vector<int> y(47);
    for(int i =0; i<47; i++){
        int suma = 0;
        for(int k =0; k<24; k++){
            if (((i-k)>=0) && ((i-k) < 24)){
                suma = suma + h[k]*x[i-k];
            }
        }
        y[i] = suma;
    }

    for(int i=0; i < y.size(); i++)
      std::cout << y.at(i) << ' ';

    ofstream outFile("x_h_conv.txt");
    if( ! outFile )	{
		cout << "Error opening file for output" << endl ;
		return -1 ;
	}
    for(int i=0; i < y.size(); i++){
        outFile << y.at(i) << endl;
    }
    outFile.close();

    return 0;
}
*/

/*
// przy uzyciu array
//399 lini asembler


#include <iostream>
#include <fstream>
#include <array>
using namespace std;

int main(){

    array<int, 47> x, h;
    ifstream f_read("x_h_sign.txt");

    x.fill(0);
    h.fill(0);

    for(int i=0; i <24; ++i){
        f_read >> x[i] >> h[i];
    }


    for(int i=0; i < x.size(); i++)
      std::cout << x.at(i) << ' ';


    array<int, 47> y;

    for(int i =0; i<47; i++){
        int suma = 0;
        for(int k =0; k<24; k++){
                suma = suma + h[k]*x[i-k];
        }
        y[i] = suma;
    }

    for(int i=0; i < y.size(); i++)
      std::cout << y.at(i) << ' ';

    ofstream outFile("x_h_conv.txt");
    if( ! outFile )	{
		cout << "Error opening file for output" << endl ;
		return -1 ;
	}
    for(int i=0; i < y.size(); i++){
        outFile << y.at(i) << endl;
    }
    outFile.close();

    return 0;
}
*/




//194 linie
#include <iostream>
#include <fstream>

using namespace std;

int main(){


    int x[47];
    int h[47];
    ifstream f_read("x_h_sign.txt");

    for(int i=0; i <47; ++i){
        x[i] = 0;
        h[i] = 0;
    }

    for(int i=0; i <24; ++i){
        f_read >> x[i] >> h[i];
    }

/*
    for(int i=0; i < x.size(); i++)
      std::cout << x.at(i) << ' ';
*/

    //list<int, 47> y;
    int y[47];

    for(int i =0; i<47; ++i){
        int suma = 0;
        for(int k =0; k<24; ++k){
            suma += h[k]*x[i-k];
        }
        y[i] = suma;
    }

    for(int i=0; i < 47; i++)
      std::cout << y[i] << ' ';

    ofstream outFile("x_h_conv.txt");
    if( ! outFile )	{
		cout << "Error opening file for output" << endl ;
		return -1 ;
	}
    for(int i=0; i < 47; i++){
        outFile << y[i] << endl;
    }
    outFile.close();

    return 0;
}
