#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <random>
#include <iterator>
#include "opencv.hpp"

using namespace cv;
using namespace std;

vector <vector<double>> cities;
vector <int> order;
int total_city_number;
random_device rd;  
mt19937 rng(rd());   

double distance_between_two(int city_1, int city_2)
{
	// the city number is the original city number, not ranged or shuffled
	double diff_x = abs(cities[city_1][0] - cities[city_2][0]);
	double diff_y = abs(cities[city_1][1] - cities[city_2][1]);
	double distance = sqrt(pow(diff_x, 2) + pow(diff_y, 2));
	return distance;
}

double calculate_total_distance()
{
	double total_distance = 0;
	for (int i = 0; i < total_city_number - 1; i++)
	{
		total_distance += distance_between_two(order[i], order[i + 1]);
	}
	total_distance += distance_between_two(order[total_city_number - 1], order[0]);
	return total_distance;
}

int get_random_number(int lower_bound, int upper_bound)
{
	// lower_bound and upper_bound are included, return int.
	uniform_int_distribution<int> uni(lower_bound, upper_bound); 
	return uni(rng);
}

void init_random_oder()
{
	for (int i = 0; i < total_city_number; i++)
	{
		order.push_back(i); // initialize a order, then shuffle it
	}
	cout << "Initial Distance :  " << calculate_total_distance() << endl;
	shuffle(begin(order), end(order), rng);
}

vector<double> read_cooridinate(string line)
{
	vector<double> result;
	istringstream iss(line);
	for (string s; iss >> s;)
	{
		result.push_back(atof(s.c_str()));
	}
	return result;
}

void nerighbours(int input, int &my_left, int &my_right)
{
	if (input != 0 && input != total_city_number - 1)
	{
		my_left = input - 1;
		my_right = input + 1;
	}
	else if (input == 0)
	{
		my_left = total_city_number - 1;
		my_right = input +1;
	}
	else if (input == total_city_number - 1)
	{
		my_left = input - 1;
		my_right = 0;
	}
}

void compute_four_distance(int the_chosen_1, int the_chosen_2,double &before, double&after)
{
	before = 0;
	after = 0;
	int chosen_1_left, chosen_1_right, chosen_2_left, chosen_2_right;
	nerighbours(the_chosen_1, chosen_1_left, chosen_1_right);
	nerighbours(the_chosen_2, chosen_2_left, chosen_2_right);

	if (abs(the_chosen_1 - the_chosen_2) == total_city_number - 1) // when its the head and tail, it is also different
	{
		if (the_chosen_1 == 0)
		{
			before += distance_between_two(order[the_chosen_1], order[chosen_1_right]);
			before += distance_between_two(order[chosen_2_left], order[the_chosen_2]);
			after += distance_between_two(order[the_chosen_2], order[chosen_1_right]);
			after += distance_between_two(order[chosen_2_left], order[the_chosen_1]);
		}
		else
		{
			before += distance_between_two(order[chosen_1_left], order[the_chosen_1]);
			before += distance_between_two(order[the_chosen_2], order[chosen_2_right]);
			after += distance_between_two(order[chosen_1_left], order[the_chosen_2]);
			after += distance_between_two(order[the_chosen_1], order[chosen_2_right]);
		}
	}
	if (abs(the_chosen_1 - the_chosen_2) == 1) //when two are neighbours, thins get tricky
	{
		if (the_chosen_1 < the_chosen_2)
		{ // when 1 is on the left side
			before += distance_between_two(order[chosen_1_left], order[the_chosen_1]);
			before += distance_between_two(order[the_chosen_2], order[chosen_2_right]);
			after += distance_between_two(order[chosen_1_left], order[the_chosen_2]);
			after += distance_between_two(order[the_chosen_1], order[chosen_2_right]);
		}
		else
		{ // when 2 is on the left side
			before += distance_between_two(order[the_chosen_1], order[chosen_1_right]);
			before += distance_between_two(order[chosen_2_left], order[the_chosen_2]);
			after += distance_between_two(order[the_chosen_2], order[chosen_1_right]);
			after += distance_between_two(order[chosen_2_left], order[the_chosen_1]);
		}
	}
	else
	{ // the common situation, have 4 neighbous to count
		before += distance_between_two(order[chosen_1_left], order[the_chosen_1]);
		before += distance_between_two(order[the_chosen_1], order[chosen_1_right]);
		before += distance_between_two(order[chosen_2_left], order[the_chosen_2]);
		before += distance_between_two(order[the_chosen_2], order[chosen_2_right]);
		after += distance_between_two(order[chosen_1_left], order[the_chosen_2]);
		after += distance_between_two(order[the_chosen_2], order[chosen_1_right]);
		after += distance_between_two(order[chosen_2_left], order[the_chosen_1]);
		after += distance_between_two(order[the_chosen_1], order[chosen_2_right]);
	}
}

void swap_two_in_order(int the_chosen_1, int the_chosen_2)
{
	// they are indexes in the order, not the order itself.
	int temp = order[the_chosen_1];
	order[the_chosen_1] = order[the_chosen_2];
	order[the_chosen_2] = temp;
}

double choose_two_and_compare(int &chosen_1, int &chosen_2)
{
	int the_chosen_1, the_chosen_2 = 0; // are index of ORDER!
	the_chosen_1 = get_random_number(0, total_city_number - 1);
	do
	{
		the_chosen_2 = get_random_number(0, total_city_number - 1);
	} while (the_chosen_2 == the_chosen_1);
	double old_distance = 0; 
	double new_distance = 0;
	// compute the old 4 distance of distance before swapping,which is 2 neighbours
	// of chosen_1 and 2 neighbours of chosen_2, and the new 4 distance. the function returns
	// the difference between them.
	compute_four_distance(the_chosen_1, the_chosen_2, old_distance, new_distance);
	chosen_1 = the_chosen_1; chosen_2 = the_chosen_2;
	return new_distance - old_distance;
}

void read_in_files(string file_string)
{
	vector <vector<double>> numbers;
	string line;
	ifstream file(file_string);
	if (file.is_open())
	{
		while (getline(file,line))
		{
			numbers.push_back(read_cooridinate(line));
		}
		file.close();
	}
	else
	{
		cout << "Cannot Open File"<<endl;
	}
	cities = numbers;
}

void see_all_orders()
{
	using iter = vector<int> ::iterator;
	cout << endl;
	for (iter it = order.begin(); it != order.end(); ++it)
		cout << *it << " ";
	cout << endl;
}


Mat drawline( int flag)
{
	Scalar gray = Scalar(200, 200, 200);
	Scalar red = Scalar(100, 100, 255);
	Mat picture(500, 500, CV_8UC3, Scalar(255, 255, 255, 0.5));
	for (int i = 0; i < total_city_number; i++)
	{
		Point center = Point(cities[i][0] * 400 + 50, cities[i][1] * 400 + 50);
		circle(picture, center, 8, Scalar(0, 0, 0), -1);
	}
	for (int i = 0; i < total_city_number-1; i++)
	{
		Point line_s = Point(cities[order[i]][0] * 400 + 50, cities[order[i]][1] * 400 + 50);
		Point line_e = Point(cities[order[i + 1]][0] * 400 + 50, cities[order[i + 1]][1] * 400 + 50);
		if (flag == 0)
			line(picture, line_s, line_e, gray, 2);
		else
			line(picture, line_s, line_e, red, 2);
	}
	Point line_s = Point(cities[order[total_city_number - 1]][0] * 400 + 50, cities[order[total_city_number - 1]][1] * 400 + 50);
	Point line_e = Point(cities[order[0]][0] * 400 + 50, cities[order[0]][1] * 400 + 50);
	if (flag == 0)
		line(picture, line_s, line_e, gray, 2);
	else
		line(picture, line_s, line_e, red, 2);
	return picture;
}
int main()
{
	read_in_files("tsp2.txt");
	total_city_number = cities.size();
	cout << "total city number:" << total_city_number << endl;
	init_random_oder();
	cout << "random initialize distance: "<<calculate_total_distance() << endl;
	cout << "initialize complete!" << endl;
	Mat picture = drawline(0);
	cv::imshow("demonstration", picture);
	cv::waitKey(1);
	//-------------------------------------------------------------------------------//
	int max_iter = 80000000
		; // nine 9s
	int iter_count = 0; 
	double completed_portion = 0;
	double temperatue = 0;

	int chosen_1, chosen_2;
	double difference = 0;
	int havent_change = 0;
	double current_distance;
	current_distance = calculate_total_distance();

	double last_accurate_distance = 0;
	double simulated_annealing_P = 0;
	double random_0_1 = 0;

	double t = cv::getTickCount();
	while (completed_portion < 1)
	{
		difference = choose_two_and_compare(chosen_1, chosen_2);
		iter_count ++ ;
		completed_portion = (double)(iter_count) / max_iter;
		temperatue = completed_portion;//log10(10 * completed_portion + 1.15);  //offset 1.1

		if (iter_count % 100 == 0)
		{
			last_accurate_distance = calculate_total_distance();
			current_distance = last_accurate_distance;
		}
		if (iter_count % 1000000 == 0)
		{
			cout << completed_portion * 100<<"% Completed" << endl;
			cout << "current distance: " << last_accurate_distance << endl;
			cout << "temperature in percent:" << temperatue * 100 << "%" << endl;
			picture = drawline(0);
			cv::imshow("demonstration", picture);
			cv::waitKey(1);
		}
		if (difference <= 0)
		{
			current_distance += difference;
			swap_two_in_order(chosen_1, chosen_2);
		}	
		else if (difference > 0)
		{
			simulated_annealing_P = 8 * (difference / current_distance) + temperatue;
			//cout << (((100)*difference) / current_distance) << endl;
			random_0_1 = ((double)get_random_number(1, 100000) / 100000);
			if (simulated_annealing_P < random_0_1)
			{
				current_distance += difference;
				swap_two_in_order(chosen_1, chosen_2);
			}
			else
			{
				havent_change++;
			}
		}
		if (havent_change > pow(total_city_number, 4.67))
			break;
	}
	cout << "Final Distance: " << calculate_total_distance() << endl;
	see_all_orders();
	picture = drawline(1);
	t = ((double)getTickCount() - t) / getTickFrequency();
	cout << "Algorithm Time" << t << "sec" << endl;
	cv::imshow("demonstration", picture);
	cv::waitKey(-1);
}
