//
//  Hill climber .cpp
//  Created by 熊之遥 on 2017/11/23.
//  Copyright © 2017年 熊之遥. All rights reserved.
//

#include <vector>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <algorithm>
#include "main.cpp"

using namespace std;
random_device rd;
mt19937 rng(rd());

class hill_climber 
{
	int iteration_limit = 0;
	int num;
	vector<CItem> items;
	double W;
	public:
		hill_climber(int iteration_limit_input, int num_input, vector<CItem> items_input, double W_input)
		{
			iteration_limit = iteration_limit_input;
			num = num_input;
			items = items_input;
			W = W_input;
		}
		CSolution * climb_one_step(CSolution* origin)
		{
			CSolution *seed = origin;
			CPopulation step_neibours(0, num, items, W, false);
			for (int i = 0; i < num;)
			{
					CSolution *one_neighbour = new CSolution(num, false);
					one_neighbour->MutateIndicator(seed, i);
					if (one_neighbour->Constrain(items, W))
					{
						i++;
						one_neighbour->Fitness(items);
						step_neibours.PutIn(one_neighbour); 
					}
			}
			step_neibours.Sort();
			CSolution * new_origin = step_neibours.PutOut();
			return new_origin;
		}
		CSolution * climb()
		{
			CSolution * init = new CSolution(num);
			for (int i = 0; i < iteration_limit; i++)
			{
				init = climb_one_step(init);
			}
			return init;
		}
};

