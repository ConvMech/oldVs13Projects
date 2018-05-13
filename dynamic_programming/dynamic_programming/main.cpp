#include<stdio.h>
#include<iostream>
#include <vector>
#include "opencv.hpp"

using namespace ::cv;
using namespace ::std;

class CItem 
{
	float value;
	int weight;
public:
	CItem(float v, float w) 
	{
		value = v;
		weight = w;
	}
	float give_value(){
		return value;
	}
	float give_weight(){
		return weight;
	}
};

struct result
{
	float best_value;
	vector<int> selection_result;
};

int max(int a, int b) { return (a > b) ? a : b; }

result knapSack(int W, vector<CItem> dictionary, int n)
{
	int i, w;
	vector<vector<result>> K(n + 1, vector<result>(W + 1));

	for (i = 0; i <= n; i++)
	{
		for (w = 0; w <= W; w++)
		{
			if (i == 0 || w == 0)
				K[i][w].best_value = 0;
			else if (dictionary[i - 1].give_weight() <= w)
			{
				float candidate = dictionary[i - 1].give_value() + K[i - 1][w - dictionary[i - 1].give_weight()].best_value;
				float without = K[i - 1][w].best_value;
				if (candidate>without)
				{
					K[i][w].best_value = candidate;
					K[i][w].selection_result = K[i - 1][w - dictionary[i - 1].give_weight()].selection_result;
					K[i][w].selection_result.push_back(i);
				}
				else{
					K[i][w].best_value = without;
					K[i][w].selection_result = K[i - 1][w].selection_result;
				}
			}
			else
				K[i][w] = K[i - 1][w];
		}
	}
	cout << K[n][W].selection_result.size() << endl;
	return K[n][W];
}

vector<CItem> initialize_dic(int val[], int wt[], int array_length)
{
	vector<CItem> dictionary;
	for (int i = 0; i < array_length; i++)
	{
		dictionary.push_back(CItem(val[i], wt[i]));
	}
	return dictionary;
}

int main()
{
	int wt[] = { 382745,
		799601,
		909247,
		729069,
		467902,
		44328,
		34610,
		698150,
		823460,
		903959,
		853665,
		551830,
		610856,
		670702,
		488960,
		951111,
		323046,
		446298,
		931161,
		31385,
		496951,
		264724,
		224916,
		169684 };
	int val[] = { 825594,
		1677009,
		1676628,
		1523970,
		943972,
		97426,
		69666,
		1296457,
		1679693,
		1902996,
		1844992,
		1049289,
		1252836,
		1319836,
		953277,
		2067538,
		675367,
		853655,
		1826027,
		65731,
		901489,
		577243,
		466257,
		369261 };
	int  W = 6404180;
	int n = sizeof(val) / sizeof(val[0]);
	cout << n << endl;
	double t = cv::getTickCount();
	
	vector<CItem> dictionary = initialize_dic(val, wt, n);
	cout << knapSack(W, dictionary, n).best_value << endl;
	vector<int> selected = knapSack(W, dictionary, n).selection_result;
	cout << selected.size() << endl;
	for (auto i = selected.begin(); i != selected.end(); ++i)
		std::cout << *i << ' ';
	
	t = ((double)getTickCount() - t) / getTickFrequency();
	cout << t << endl;
	system("pause");
	return 0;
}