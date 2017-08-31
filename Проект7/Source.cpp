#define _USE_MATH_DEFINES
#include <iostream>
#include <time.h>
#include <vector>
#include <cmath>
using namespace std;
double BinDec(vector<char> bin, int size)
{
	int rez = 0, i;
	for (int i = size - 1; i >= 0; --i)
		rez += (bin[i] - 48) << i;
	return rez;
}
vector<char> DecBin(int chislo, int size)
{
	vector<char> bin(size);
	for (int i = size - 1; i >= 0; --i)
		bin[i] = ((chislo >> i) & 1) + 48;
		
	return bin;
}
double Heming(vector<char> Obj1, vector<char> Obj2)
{
	int ans = 0;
	for (int i = 0; i < Obj1.size(); i++) {
		if (Obj1[i] != Obj2[i])
			ans++;
	}
	return ans;
}
void Cross(vector<vector<char>> Obj_Parents, int N, vector<vector<double>>Sign)
{
	vector<int> key(Obj_Parents.size());
	vector<char> temp(N);
	vector<vector<char>> Obj_Kids;
	vector<char> Kids1(N); vector<char> Kids2(N);
	
	int A, KeyI, KeyJ;
	for (int i = 0; i < Obj_Parents.size(); i++) {
		if (key[i] == 1) continue;
		temp = Obj_Parents[i]; key[i] = 1; KeyI = i;
		vector<int> amount(Obj_Parents.size());
		for (int j = 0; j < Obj_Parents.size(); j++) {
			if (key[j] == 1) {
				amount[j] = 0;
			}
			else {
				amount[j] = Heming(temp, Obj_Parents[j]);
			}		
		}
		A = -10000;
		for (int j = 0; j < amount.size(); j++) {
			if (A < amount[j]) {
				A = amount[j];
				KeyJ = j;
			}
		}
		key[KeyJ] = 1;
		for (int j = 0; j < N / 2; j++) {
			Kids1[j] = Obj_Parents[KeyI][j];
			Kids2[j] = Obj_Parents[KeyJ][j];
		}
		for (int j = N / 2; j < N; j++) {
			Kids1[j] = Obj_Parents[KeyJ][j];
			Kids2[j] = Obj_Parents[KeyI][j];
		}
		Obj_Kids.push_back(Kids1);
		Obj_Kids.push_back(Kids2);
	}
	//return Obj_Kids;
}

vector<vector<char>> Mutation(vector<vector<char>> Obj_Parents, int N)
{
	for (int i = 1; i < Obj_Parents.size(); i++) {
		for (int j = 0; j < N; j++) {
			if (rand() % 2 == 1)
				Obj_Parents[i][j] = '1' - Obj_Parents[i][j] + 48;
		}
	}
	return Obj_Parents;
}

double FF(double x, double y)
{	
	//return (double)(x*x);
	//return (double)(x*x - 10 * cos(2 * M_PI* x));
	//return abs(x);
	//return 2 * x * x - (double)(1.05 * x*x*x*x) + (double)((x*x*x*x*x*x) / 6) + x*y + y*y;
	//return 2 * x*x - 1.05*x*x*x*x + (x*x*x*x*x*x) / 6 + x*y + y*y;
	return  -(y + 47) * sin(sqrt(abs(x / 2 + (y + 47)))) - x * sin(sqrt(abs(x - (y + 47))));
	//return -cos(x)*cos(y)*exp(-((x - M_PI)*(x - M_PI) + (y - M_PI)*(y - M_PI)));
}
int main()
{
	setlocale(LC_ALL, "Russian");
	int Nobj, Ngen, N = 14, Max_Iter = 200;
	int a = -100, b  = 100;
	int c = pow(2, N) / (b - a);
	cin >> Nobj >> Ngen;
	vector<vector<double>> Obj(Nobj, vector<double>(Ngen, 0));
	vector<vector<char>> Obj_Parents(Nobj, vector<char>(N*Ngen, 0));
	vector<vector<char>> Obj_Kids(Nobj, vector<char>(N*Ngen, 0));
	vector<vector<char>> Kids_Parents(Nobj*2, vector<char>(N*Ngen, 0));
	vector<vector<double>> Numerical(Nobj * 2, vector<double>(Ngen, 0));
	vector<vector<double>> Sign(Nobj * 2, vector<double>(Ngen, 0));
	vector<double> value_ff(Nobj * 2);
	vector<vector<char>> Obj_Ans(Nobj, vector<char>(N*Ngen, 0));
	vector<vector<double>> Obj_Ans_Num(Nobj, vector<double>(Ngen, 0));
	double exit = 10000000;
	srand(time(NULL));
	for (int i = 0; i < Nobj; i++) {
		for (int j = 0; j < Ngen; j++) {
			int PP = rand() % (b*1000);
			double PP1 =  (double)(PP) / 1000;
			double M = a  + rand() % ((b - a));
			double Z = abs(a  + rand() % ((b - a)));
			if (Z == 0) {
				Z = M; 
			}
				Obj[i][j] = (double) (M / Z) * 10000;
			//cout << Obj[i][j] << " ";
		}
		//cout << endl;
	}
	vector<char> temp(N);
	int Temp;
	for (int i = 0; i < Nobj; i++) {
		Temp = 0;
		for (int j = 0; j < Ngen; j++) {
			temp = DecBin(Obj[i][j], N);
			for (int k = 0; k < N; k++) {
				Obj_Parents[i][k + Temp] = temp[k];
			}
			Temp += N;
		}
	}
	//int coll = 0, rav = 0;
	for (int i = 0; i < Max_Iter; i++)
	{
		//Скрещивание
		 Cross(Obj_Parents, N*Ngen, Sign); 
		//Объединение детей и родителей
		for (int j = 0; j < Nobj; j++) { 
			for (int k = 0; k < N*Ngen; k++) {
				Kids_Parents[j][k] = Obj_Parents[j][k];
				Kids_Parents[j + Nobj][k] = Obj_Kids[j][k];
			}
		}
		//cout << endl;
		vector<vector<char>> Temp1(Ngen, vector<char>(N, 0));
			//Перевод из бинарной в алгебраическую структуру
			for (int j = 0; j < Nobj * 2; j++) {
				for (int k = 0; k < Ngen; k++) {
					for (int l = 0 ; l < N; l++) {
						Temp1[k][l] = Kids_Parents[j][l+N*k];
					}
				}
				for (int k = 0; k < Ngen; k++) {
					double C = (double)(BinDec(Temp1[k],N));
					Numerical[j][k] = C / 10000;
					//cout << BinDec(Temp1[k], N) / 10000 << " ";
				}
				//cout << endl;
			}

			/*for (int j = 0; j < Nobj * 2; j++) {
				for (int k = 0; k < Ngen; k++) {
					if (rand() % 2 == 0) Numerical[j][k] *= -1;
				}
			}*/

			for (int j = 0; j < Nobj ; j++) {
				value_ff[j] = 0;
			}
			int PP = 20;
			//Просчет FitnesFunction
			for (int j = 0; j < Nobj * 2; j++) {
					value_ff[j] = FF(Numerical[j][0], Numerical[j][1]);
			}
			// -100 100  = -100 + (rand() % (200*10000)) / 10000;
			//Сортировка
			for (int j = 0; j < Nobj * 2; j++) {
				for (int k = j + 1; k < Nobj*2; k++) {
					if (value_ff[j] > value_ff[k]) {
						double temp = value_ff[k];
						value_ff[k] = value_ff[j];
						value_ff[j] = temp;
						vector<char> temp1(N*Ngen);
						temp1 = Kids_Parents[k];
						Kids_Parents[k] = Kids_Parents[j];
						Kids_Parents[j] = temp1;
					}
				}
			}

			for (int j = 0; j < Nobj; j++) {
				for (int k = 0; k < N*Ngen; k++) {
					Obj_Ans[j][k] = Kids_Parents[j][k];
				}
			}
			// Мутация
			Obj_Parents = Mutation(Obj_Ans, N*Ngen);
			Obj_Parents = Mutation(Obj_Parents, N*Ngen);
			//Obj_Parents = Mutation(Obj_Parents, N*Ngen);
			vector<vector<char>> Temp2(Ngen, vector<char>(N, 0));
			for (int j = 0; j < Nobj; j++) {
				for (int k = 0; k < Ngen; k++) {
					for (int l = 0; l < N; l++) {
						Temp2[k][l] = Kids_Parents[j][l + N*k];
					}
				}
				for (int k = 0; k < Ngen; k++) {
					Obj_Ans_Num[j][k] = BinDec(Temp2[k], N) / 10000;
				}
			}
				
			for (int j = 0; j < Nobj; j++) {
				for (int k = 0; k < Ngen; k++) {
					cout << Obj_Ans_Num[j][k] << " ";
				}
				cout << value_ff[j] << " ";
				cout << endl;
			}
			cout << endl;
	}
	system("pause");
	return 0;
}