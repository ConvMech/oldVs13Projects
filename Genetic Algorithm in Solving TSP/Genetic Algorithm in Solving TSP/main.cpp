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

vector<pair<int, float>> performance;
vector <vector<double>> cities;   // all the cities's location information
vector <int> order;  // the order of the city 
int total_city_number;
random_device rd;
mt19937 rng(rd());

void write_performance( string file, float final_distance, float t)
{
	ofstream myfile;
	myfile.open(file);
	for (int i = 0; i < performance.size(); i++)
	{
		myfile << performance[i].first << "-" << performance[i].second << "\n";
	}
	for (int i = 0; i < order.size(); i++)
	{
		myfile << order[i]<<" ";
	}
	myfile << "\n" << "final_distance" << final_distance << endl;
	myfile << "\n" << "total_time" << t << endl;
	myfile.close();
}

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
	order.clear();
	for (int i = 0; i < total_city_number; i++)
	{
		order.push_back(i); // initialize a order, then shuffle it
	}
	//cout << "Initial Distance :  " << calculate_total_distance() << endl;
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
		my_right = input + 1;
	}
	else if (input == total_city_number - 1)
	{
		my_left = input - 1;
		my_right = 0;
	}
}

void compute_four_distance(int the_chosen_1, int the_chosen_2, double &before, double&after)
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
		while (getline(file, line))
		{
			numbers.push_back(read_cooridinate(line));
		}
		file.close();
	}
	else
	{
		cout << "Cannot Open File" << endl;
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


Mat drawline(int flag)
{
	Scalar gray = Scalar(200, 200, 200);
	Scalar red = Scalar(100, 100, 255);
	Mat picture(1000, 1000, CV_8UC3, Scalar(255, 255, 255, 0.5));
	for (int i = 0; i < total_city_number; i++)
	{
		Point center = Point(cities[i][0] * 300 + 500, cities[i][1] * 300 + 500);
		circle(picture, center, 8, Scalar(0, 0, 0), -1);
	}
	for (int i = 0; i < total_city_number - 1; i++)
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

void simulated_annealing( )
{
	int max_iter = 8000000
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

	while (completed_portion < 1)
	{
		difference = choose_two_and_compare(chosen_1, chosen_2);
		iter_count++;
		completed_portion = (double)(iter_count) / max_iter;
		temperatue = completed_portion;//log10(10 * completed_portion + 1.15);  //offset 1.1

		if (iter_count % 100 == 0)
		{
			last_accurate_distance = calculate_total_distance();
			current_distance = last_accurate_distance;
		}
		if (iter_count % 90000 == 0)
		{
			cout << completed_portion * 100 << "% Completed" << endl;
			cout << "current distance: " << last_accurate_distance << endl;
			cout << "temperature in percent:" << temperatue * 100 << "%" << endl;
			performance.push_back(make_pair(iter_count, last_accurate_distance));
			Mat picture = drawline(0);
			cv::imshow("demonstration", picture);
			cv::waitKey(1);
		}
		if (difference >= 0) // <= for minimum
		{
			current_distance += difference;
			swap_two_in_order(chosen_1, chosen_2);
		}
		else if (difference < 0)// >for minimum
		{
			// no (-1)*for minimum
			simulated_annealing_P = 8 * ((-1)*difference / current_distance) + temperatue;
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
}
vector<int> mutation(double rate, vector<int> one_chromosome )
{
	double random_0_1 = ((double)get_random_number(1, 100000) / 100000);
	if (random_0_1 < rate)
	{
		int the_chosen_1, the_chosen_2 = 0;
		the_chosen_1 = get_random_number(0, total_city_number - 1);
		do
		{
			the_chosen_2 = get_random_number(0, total_city_number - 1);
		} while (the_chosen_2 == the_chosen_1);
		int temp = one_chromosome[the_chosen_1];
		one_chromosome[the_chosen_1] = one_chromosome[the_chosen_2];
		one_chromosome[the_chosen_2] = temp;
	}
	return one_chromosome;
}

int tourment_selection(int k, vector<pair<float, int>> fitness_and_index)
{
	int best_index = -1;
	for (int i = 0; i < k; i++)
	{
		int random = get_random_number(0, fitness_and_index.size() - 1);
		if (best_index == -1 || fitness_and_index[random].first < fitness_and_index[best_index].first)
			best_index = random;
	}
	return best_index;
}

vector< vector<int>> cross_breeding(vector<int>chomosome1, vector<int>chomosome2)
{
	int the_chosen_1, the_chosen_2 = 0; // get_the_position_of_the_cut
	vector< vector<int>> two_childern;
	vector<int> child_1(total_city_number, 0);
	vector<int>child_2(total_city_number, 0);
	// get_two_random_numbers;
	the_chosen_1 = get_random_number(0, total_city_number - 1);
	do
	{
		the_chosen_2 = get_random_number(0, total_city_number - 1);
	} while (the_chosen_2 == the_chosen_1);
	int smaller, bigger = 0;
	if (the_chosen_1 < the_chosen_2)	{
		smaller = the_chosen_1; bigger = the_chosen_2;
	}
	else{
		smaller = the_chosen_2; bigger = the_chosen_1;
	}
	//smaller = 2; bigger = 3;
	//cout << "the_chosen_positions_are" << smaller << "and" << bigger << endl;
	vector<int> chromo1_hash(total_city_number, 0); // hash table to write the OX confilict
	vector<int> chromo2_hash(total_city_number, 0);
	vector<int> chromo1_remain, chromo2_remain;
	for (int i = smaller; i <= bigger; i++)
	{
		child_1[i] = chomosome1[i];
		chromo1_hash[chomosome1[i]] = 1;
		child_2[i] = chomosome2[i];
		chromo2_hash[chomosome2[i]] = 1;
	}
	for (int i = bigger + 1; i <= total_city_number - 1; i++) // add the right remaining part to remian
	{
		chromo1_remain.push_back(chomosome1[i]);
		chromo2_remain.push_back(chomosome2[i]);
	}
	for (int i = 0; i <= smaller - 1; i++) // now the left remaining part
	{
		chromo1_remain.push_back(chomosome1[i]);
		chromo2_remain.push_back(chomosome2[i]);
	}
	// now the central not-duplicate part
	for (int i = smaller; i <= bigger; i++)
	{
		if (chromo2_hash[chomosome1[i]] == 0)
			// meaning that this number in the center does not exit in chromo_2
		{
			chromo1_remain.push_back(chomosome1[i]);
		}
		if (chromo1_hash[chomosome2[i]] == 0)
			// same thing
		{
			chromo2_remain.push_back(chomosome2[i]);
		}
	}
	//now_let's generate the children!
	int j_1 = 0,j_2 = 0;
	for (int i = bigger + 1; i <= total_city_number - 1; i++) // add the right remaining part to remian
	{
		if (chromo1_hash[chromo2_remain[j_1]] == 0)
		{
			child_1[i] = chromo2_remain[j_1];
			j_1++;
		}
		else
		{
			while (chromo1_hash[chromo2_remain[j_1]] != 0)
				j_1++;
			child_1[i] = chromo2_remain[j_1];
			j_1++;
		}
		if (chromo2_hash[chromo1_remain[j_2]] == 0)
		{
			child_2[i] = chromo1_remain[j_2];
			j_2++;
		}
		else
		{
			while (chromo2_hash[chromo1_remain[j_2]] != 0)
				j_2++;
			child_2[i] = chromo1_remain[j_2];
			j_2++;
		}
	}
	for (int i = 0; i <= smaller - 1; i++) // add the left remaining part to remian
	{
		if (chromo1_hash[chromo2_remain[j_1]] == 0)
		{
			child_1[i] = chromo2_remain[j_1];
			j_1++;
		}
		else
		{
			while (chromo1_hash[chromo2_remain[j_1]] != 0)
				j_1++;
			child_1[i] = chromo2_remain[j_1];
			j_1++;
		}
		if (chromo2_hash[chromo1_remain[j_2]] == 0)
		{
			child_2[i] = chromo1_remain[j_2];
			j_2++;
		}
		else
		{
			while (chromo2_hash[chromo1_remain[j_2]] != 0)
				j_2++;
			child_2[i] = chromo1_remain[j_2];
			j_2++;
		}
	}
	two_childern.push_back(child_1); two_childern.push_back(child_2);
	return two_childern;
}

float calculate_fitness(vector<int> chromosome)
{
	float total_distance = 0;
	for (int i = 0; i < total_city_number - 1; i++)
	{
		total_distance += distance_between_two(chromosome[i], chromosome[i + 1]);
	}
	total_distance += distance_between_two(chromosome[total_city_number - 1], chromosome[0]);
	return total_distance;
}

vector< vector<int>>breed_mutate_selection(int eltism_num, int crossover_num, double mutation_rate,  int tourment_value, vector< vector<int>> chromosomeset)
{
	vector< vector<int>> new_generation;
	vector<pair<float, int>> fitness_and_index; // fitness_on_the_first,index_on_the_second
	for (int i = 0; i < chromosomeset.size(); i++)
	{
		float this_fittness = calculate_fitness(chromosomeset[i]);
		fitness_and_index.push_back(make_pair(this_fittness, i));
		//cout << "number" << i<< " fitness" << this_fittness << endl;
	}
	sort(fitness_and_index.begin(), fitness_and_index.end(), [](const std::pair<float, int> &left, const std::pair<float, int> &right)
	{
		return left.first > right.first;   // for shortest path < 
	});

	int new_g_i = 0;
	vector< vector<int>> breeding_parents;
	vector< vector<int>> elites;
	while (new_g_i < eltism_num) //preserve the elite of the chromoset
	{
		vector<int> elite_itself = chromosomeset[fitness_and_index[new_g_i].second];
		elites.push_back(mutation(mutation_rate, elite_itself));
		new_g_i++;
	}
	vector<int > selected_parents; //breed a number of parents
	for (int j = 0; j < eltism_num + crossover_num; j++) // try to use selection function next
	{
		selected_parents.push_back(j); // take_these_parents
	}
	shuffle(begin(selected_parents), end(selected_parents), rng); // now_we_have_shuffled_parents
	for (int i = 0; i < selected_parents.size()/2; i++)
	{
		new_g_i= new_g_i + 2;
		vector< vector<int>> two_children = cross_breeding(chromosomeset[fitness_and_index[selected_parents[i]].second],
			chromosomeset[fitness_and_index[selected_parents[i + 1]].second]);
		new_generation.push_back(mutation(mutation_rate, two_children[0]));
		new_generation.push_back(mutation(mutation_rate, two_children[1]));
	}
	while (new_g_i < chromosomeset.size()) // the_remaining_chroms
	{
		vector<int> looser = chromosomeset[fitness_and_index[new_g_i].second];
		new_generation.push_back(mutation(mutation_rate, looser));
		new_g_i++;
	}

	for (int i = 0; i < new_generation.size(); i++)
	{
		float this_fittness = calculate_fitness(new_generation[i]);
		fitness_and_index.push_back(make_pair(this_fittness, chromosomeset.size() + i));
		//cout << "number" << i<< " fitness" << this_fittness << endl;
	}
	sort(fitness_and_index.begin(), fitness_and_index.end(), [](const std::pair<float, int> &left, const std::pair<float, int> &right)
	{
		return left.first < right.first;
	});

	vector< vector<int>> next_generation;
	vector<int>good_chromosome;
	int i = 0;
	for (int i = 0; i < elites.size(); i++)
		next_generation.push_back(elites[i]);
	while (next_generation.size() < chromosomeset.size())
	{
		// tourment_selection
		int index = tourment_selection(tourment_value, fitness_and_index);
		if (index < chromosomeset.size())
			good_chromosome = chromosomeset[index];
		else
			good_chromosome = new_generation[index - chromosomeset.size()];
		i++;
		next_generation.push_back(good_chromosome);
	}
	
	return next_generation;
}

vector< vector<int>> generate_initial_c_set( int P)
{
	vector< vector<int>> chomosomeset;
	for (int i = 0; i < P; i++)
	{
		vector<int > chomosomeset_sub;
		for (int j = 0; j < total_city_number; j++)
		{
			chomosomeset_sub.push_back(j); // initialize a order, then shuffle it
		}
		shuffle(begin(chomosomeset_sub), end(chomosomeset_sub), rng);
		chomosomeset.push_back(chomosomeset_sub);
	}
	return chomosomeset;
}

void see_all_chromosme(int P, vector< vector<int>> chromosomeset)
{
	vector<int > chromosomeset_sub;
	for (int i = 0; i < P; i++)
	{
		chromosomeset_sub = chromosomeset[i];
		using iter = vector<int> ::iterator;
		cout <<"CHROMOSOME NUMBER: "<<i<< endl;
		for (iter it = chromosomeset_sub.begin(); it != chromosomeset_sub.end(); ++it)
			cout << *it << " ";
		cout << endl;
	}
}

void genetic_algorithm_order_crossover() // P would be the number of GA each revolution
{
	int P = 36;
	int loop = 1000000;
	vector< vector<int>> chromosomeset;
	chromosomeset = generate_initial_c_set(P);
	/*see_all_chromosme(P, chromosomeset);
	chromosomeset = cross_breeding(chromosomeset[0], chromosomeset[1]);
	see_all_chromosme(P, chromosomeset);*/
	int real_loop = (1000000 - loop);
	while (real_loop < 8000000)
	{
		chromosomeset = breed_mutate_selection(2, 34, 0.02,2,chromosomeset);
		loop--;
		real_loop = (1000000 - loop)*P;
		if (real_loop % 10000 == 0)
		{
			cout << real_loop;
			float distance_now = calculate_fitness(chromosomeset[0]);
			cout << "  " << distance_now << endl;
			performance.push_back(make_pair(real_loop, distance_now));
			order = chromosomeset[0];
			Mat picture = drawline(0);
			cv::imshow("demonstration", picture);
			cv::waitKey(1);
			//see_all_chromosme(P, chromosomeset);
		}
	}
	
}

void random_algorithm()
{
	int loop = 0;
	float distance = -1;
	vector <int> good_order;
	while (loop < 8000000)
	{
		loop++;
		init_random_oder();
		float new_distance = calculate_total_distance();
		if (new_distance > distance || distance == -1)  // < fpr minimum
		{
			distance = new_distance;
			good_order = order;
		}
		if (loop % 90000 == 0)
		{
			cout << "loop:" << loop << "  " << distance << endl;
			performance.push_back(make_pair(loop, distance));
			Mat picture = drawline(0);
			cv::imshow("demonstration", picture);
			cv::waitKey(1);
		}
	}
	order = good_order;
}

int main()
{
	read_in_files(" .txt");
	total_city_number = cities.size();
	cout << "total city number:" << total_city_number << endl;
	init_random_oder();
	cout << "random initialize distance: " << calculate_total_distance() << endl;
	cout << "initialize complete!" << endl;
	performance.push_back(make_pair(0, calculate_total_distance()));
	Mat picture = drawline(0);
	cv::namedWindow("demonstration");
	cv::imshow("demonstration", picture);
	cv::waitKey(1);
	double t = cv::getTickCount();
	
	//-------------------------------------------------------------------------------//
	//random_algorithm();
	simulated_annealing();
	//genetic_algorithm_order_crossover();
	//-------------------------------------------------------------------------------//
	cout << "Final Distance: " << calculate_total_distance() << endl;
	see_all_orders();
	picture = drawline(1);
	t = ((double)getTickCount() - t) / getTickFrequency();
	write_performance("result_store.txt", calculate_total_distance(),t);
	cout << "Algorithm Time" << t << "sec" << endl;
	cv::imshow("demonstration", picture);
	cv::waitKey(-1);
}