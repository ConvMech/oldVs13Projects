#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "opencv.hpp"
#include <cmath>

using namespace ::cv;
using namespace ::std;

random_device rd;
mt19937 rng(rd());

vector <vector<double>> test_xy;
vector<pair<int, float>> performance;

int record_frequncy = 10;
int maxim_layer = 4;
float penalty = 0.02; //0.02
float random_l = -50;   // the lower range to random a number
float random_r = 50;    // the upper range to random a number

enum type { plus, minus, mul, divid, number, symbol};
type type_array[7] = { type::plus, type::minus, type::mul, type::divid, type::number, type::symbol, type::symbol };


struct tree_node
{
	int number;  
	char symbol;   
	type gender;
	struct tree_node * leftchild;
	struct tree_node * rightchild;
	int depth;
};
pair<tree_node*, float> best;
tree_node *test_tree = new tree_node();
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

void write_performance(string file, float final_distance, float t)
{
	ofstream myfile;
	myfile.open(file);
	for (int i = 0; i < performance.size(); i++)
	{
		myfile << performance[i].first << "-" << performance[i].second << "\n";
	}
	myfile << "\n" << "final_distance" << final_distance << endl;
	myfile << "\n" << "total_time" << t << endl;
	myfile.close();
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
	test_xy = numbers;
}

void print_operator(tree_node * &T)
{
	type gender = T->gender;
	if (gender == type::plus)
		cout << '+' ;
	else if (gender == type::minus)
		cout << '-' ;
	else if (gender == type::mul)
		cout << '*' ;
	else if (gender == type::divid)
		cout << '/' ;
	else if (gender == type::symbol)
		cout << T->symbol;
	else if (gender == type::number)
		cout << "a_number" << endl;
}
int random_int(int lower_bound, int upper_bound)
{
	uniform_int_distribution<int> uni(lower_bound, upper_bound);
	return uni(rng);
}

type random_type( int flag = 0)
{
	if (flag == 0)
	{
		uniform_int_distribution<int> uni(0, 70);
		return type_array[(int)(uni(rng) / 10)];
	}
	else
	{
		uniform_int_distribution<int> uni(0, 40);
		return type_array[(int)(uni(rng) / 10)];
	}
}

float random_float(float lower, float upper)
{
	uniform_int_distribution<int> uni(-100, 100);
	float kernel = (float)(uni(rng) / (float)200);
	return ((upper - lower) / 2 + (upper - lower)*kernel);
}

int random_boolean()
{
	uniform_int_distribution<int> uni(0, 200);
	return type_array[(int)(uni(rng) / 100)];
}

int random_prob(float probablity)
{
	int accuracy = 1000000;
	uniform_int_distribution<int> uni(0, accuracy);
	int random = (int)(uni(rng));
	if ((probablity*accuracy) > random)
		return 1;
	else
		return 0;
}

void grow_a_tree( tree_node * &T,int depth)
{
	T = new tree_node();
	T->depth = depth;
	if (T->depth == maxim_layer)// reached limit
	{
		if (random_boolean())
		{
			T->gender = type::number;
			T->number = random_float(-1, 1);
			return;
		}
		else
		{
			T->gender = type::symbol;
			T->symbol = 'x';
			//if (random_boolean())
			//	T->symbol= 'x';
			//else
			//	T->symbol = 'y';
			return;
		}
	}
	
	if (T->depth == 0) // that is a root
	{
		T->gender = random_type(1);  // assign a operator
		while (T->gender == type::number || T->gender == type::symbol)
		{
			T->gender = random_type(1);
		}
	}
	else
	{
		T->gender = random_type();
		if (T->gender == type::number)
		{
			//T->number = random_float(random_l, random_r);
			T->number = random_int(random_l, random_r);
			return;
		}
		else if (T->gender == type::symbol)
		{
			T->symbol = 'x';
			//if (random_boolean())
			//	T->symbol= 'x';
			//else
			//	T->symbol = 'y';
			return;
		}
	}
	grow_a_tree(T->leftchild, depth + 1);
	grow_a_tree(T->rightchild, depth + 1);
}

void print_tree( tree_node * &T, string flag = "print")
{
	if (T) 
	{
		if (flag == "print"&&T->leftchild != NULL&&T->depth != 0)
			cout << "(";
		print_tree(T->leftchild, flag);

		if (T->gender != type::number)
		{
			//cout << " d" << T->depth<<": ";
			print_operator(T);
			//cout << T->depth;
		}
		else
		{ 
			//cout << " d" << T->depth << ": ";
			cout << T->number;
			//cout << T->depth;
		}
		print_tree(T->rightchild, flag);
		if (flag == "print"&&T->rightchild != NULL&&T->depth != 0)
			cout << ")";
	}
	else cout << "";
}

float simple_cal(float num1, float num2, type gender)
{
	if (gender == type::plus)
		return (num1 + num2);
	else if (gender == type::minus)
		return (num1 - num2);
	else if (gender == type::mul)
		return (num1 * num2);
	else if (gender == type::divid)
	{
		if (num2 != 0)
			return (num1 / num2);
		else
			return 1;
	}
}

float compute_tree(tree_node * &T, float real_value_x, float real_value_y) // calculate value
{
	float result;
	if (T)
	{
		if (T->gender == type::symbol)
		{
			if (T->symbol == 'x')
				return real_value_x;
			else if (T->symbol == 'y')
				return real_value_y;
		}
		else if (T->gender == type::number)
		{
			return T->number;
		}
		else
		{
			float left = compute_tree(T->leftchild, real_value_x, real_value_y);
			float right = compute_tree(T->rightchild, real_value_x, real_value_y);
			result = simple_cal(left, right, T->gender);
			return result;
		}
	}
}

int update_depth(tree_node * &T, int now_depth) // init a depth cout start
{
	if (T == NULL){
		return (now_depth - 1);
	}
	T->depth = now_depth;
	int nLeft = update_depth(T->leftchild, now_depth + 1);
	int nRight = update_depth(T->rightchild, now_depth + 1);
	if (nLeft > nRight)
		return nLeft;
	else
		return nRight;
}

// to gain the equal chance, the now number should be initiated as 1 
int random_choose_node(tree_node * &T,int now_number, tree_node * &selected)
{
	if (T == NULL){
		return (now_number - 1);  // now number is the flag to count the number of traversed nodes
	}
	int random_kernel = random_int(1, now_number);
	if (random_kernel == 1)
		if (T->gender != type::symbol && T->gender != type::number)
			selected = T;

	now_number = random_choose_node(T->leftchild, now_number + 1, selected);
	now_number = random_choose_node(T->rightchild, now_number + 1, selected);
	return now_number;
}

int random_choose_number(tree_node * &T, int now_number, tree_node * &selected)
{
	if (T == NULL){
		return (now_number);  // now number is the flag to count the number of traversed nodes
	}
	if (T->gender == type::number)
	{
		now_number++;
		int random_kernel = random_int(1, now_number);
		if (random_kernel == 1)
			selected = T;
	}
	now_number = random_choose_number(T->leftchild, now_number, selected);
	now_number = random_choose_number(T->rightchild, now_number, selected);
	return now_number;
}

tree_node* clone_a_tree(tree_node * &T)
{
	if (T != NULL)
	{
		tree_node *clone_node = new tree_node();
		clone_node->depth = T->depth;
		clone_node->gender= T->gender;
		clone_node->number = T->number;
		clone_node->symbol = T->symbol;
		clone_node->leftchild = clone_a_tree(T->leftchild);
		clone_node->rightchild = clone_a_tree(T->rightchild);
		return clone_node;
	}
	else
		return NULL;
}

void delete_tree(tree_node * &t)
{
	if (t == NULL)
	{
		return;
	}
	delete_tree(t->leftchild);
	delete_tree(t->rightchild);
	delete t;
}

void mutation(tree_node * &T1)
{
	tree_node * t1selected;
	int t1_number = random_choose_node(T1, 1, t1selected); // get the chosen node
	int current_depth = t1selected->depth;
	tree_node * new_sub_tree = new tree_node(); // create a node to grow a new tree
	grow_a_tree(new_sub_tree, current_depth + 1);
	if (random_boolean())
	{
		// left child delete, replace
		delete_tree(t1selected->leftchild);
		t1selected->leftchild = new_sub_tree;
	}
	else
	{
		delete_tree(t1selected->rightchild);
		t1selected->rightchild = new_sub_tree;
	}
	if (!random_prob(1.0 / t1_number))
	{
		T1->gender = random_type(1);
		while (T1->gender == type::number || T1->gender == type::symbol)
		{
			T1->gender = random_type(1);
		}
	}
}

void mutation_number(tree_node * &Tn)
{
	tree_node * num_selected;
	int t1_number = random_choose_number(Tn, 1, num_selected);
	if (num_selected == NULL)
		return;
	if (num_selected->gender != type::number)
		return;
	else
	{
		if (random_boolean())
			num_selected->number = num_selected->number + 1;
		else
			num_selected->number = num_selected->number - 1;

	}
}

vector<tree_node*>cross_over(tree_node * &T1, tree_node * &T2)
{
	tree_node * t1selected; tree_node * t2selected;
	tree_node * t1cloned = clone_a_tree(T1);
	tree_node * t2cloned = clone_a_tree(T2);

	int t1_number = random_choose_node(t1cloned, 1, t1selected);
	int t2_number = random_choose_node(t2cloned, 1, t2selected);
	if (random_boolean()) // t1choose_left
	{
		tree_node * t1_subtree = clone_a_tree(t1selected->leftchild);
		delete_tree(t1selected->leftchild);
		if (random_boolean()) // t2choose_left
		{
			tree_node * t2_subtree = clone_a_tree(t2selected->leftchild);
			delete_tree(t2selected->leftchild);
			t2selected->leftchild = t1_subtree;
			t1selected->leftchild = t2_subtree;
		}
		else //t2choose_right
		{
			tree_node * t2_subtree = clone_a_tree(t2selected->rightchild);
			delete_tree(t2selected->rightchild);
			t2selected->rightchild = t1_subtree;
			t1selected->leftchild = t2_subtree;
		}
	}
	else
	{
		tree_node * t1_subtree = clone_a_tree(t1selected->rightchild);
		delete_tree(t1selected->rightchild);
		if (random_boolean()) // t2choose_left
		{
			tree_node * t2_subtree = clone_a_tree(t2selected->leftchild);
			delete_tree(t2selected->leftchild);
			t2selected->leftchild = t1_subtree;
			t1selected->rightchild = t2_subtree;
		}
		else //t2choose_right
		{
			tree_node * t2_subtree = clone_a_tree(t2selected->rightchild);
			delete_tree(t2selected->rightchild);
			t2selected->rightchild = t1_subtree;
			t1selected->rightchild = t2_subtree;
		}
	}
	vector<tree_node*> children;
	children.push_back(t1cloned); children.push_back(t2cloned);
	return children;
}

vector<tree_node*>cross_over_with_limit(tree_node * &T1, tree_node * &T2)
{
	vector<tree_node*> children = cross_over(T1, T2);
	int new_t1_depth = update_depth(children[0], 0);
	int new_t2_depth = update_depth(children[1], 0);
	while (new_t1_depth > maxim_layer || new_t2_depth > maxim_layer)
	{
		delete_tree(children[0]);
		delete_tree(children[1]);
		children = cross_over(T1, T2);
		new_t1_depth = update_depth(children[0], 0);
		new_t2_depth = update_depth(children[1], 0);
	}
	return children;
}

int count_a_tree(tree_node * &T,int num = 0)
{
	if (T == NULL)
		return num;
	else
	{
		int left = count_a_tree(T->leftchild, num);
		int right = count_a_tree(T->rightchild, num);
		return left + right + 1;
	}
}

float get_fitness_value(tree_node * &T,int p = 1)
{
	float total_fitness_value = 0;
	int tri_penalty = 0;
	if (compute_tree(T, 10000000000, 10000000000) < 1)
	{
		tri_penalty = 200;
	}
	if (abs(compute_tree(T, 10, 10) - compute_tree(T, 10, 10000000000)) == 0)
	{
		tri_penalty = 100;
	}
	if (abs(compute_tree(T, 10000000000, 10) - compute_tree(T, 10, 10)) == 0)
	{
		tri_penalty = 100;
	}
	for (int i = 0; i < test_xy.size(); i++)
	{
		float x = test_xy[i][0]; float y = test_xy[i][1];
		float equation = compute_tree(T, x, y);
		total_fitness_value += abs(y - equation);
	}
	if (total_fitness_value == 0)
		tri_penalty = 100;
	float weighted_num = penalty * count_a_tree(T);
	if (p == 0)
		weighted_num = 0;
	tri_penalty = 0; // cases when use only x
	return total_fitness_value + weighted_num + tri_penalty;
}

vector<pair<tree_node*, float>> init_gp(int pool)
{
	vector<pair<tree_node*,float>> first_generation;
	for (int i = 0; i < pool; i++)
	{
		tree_node * rootNode = new tree_node();
		grow_a_tree(rootNode, 0);
		float fitness = get_fitness_value(rootNode);
		first_generation.push_back(make_pair(rootNode, fitness));
	}
	return first_generation;
}

int tourment_selection(int k, vector< pair<tree_node*, float>> caled_selection_pool)
{
	int best_index = -1; // init_best_index
	for (int i = 0; i < k; i++)
	{
		int random = random_int(0, caled_selection_pool.size() - 1);
		if (best_index == -1 || caled_selection_pool[random].second < caled_selection_pool[best_index].second)
			best_index = random;
	}
	return best_index;
}

Mat drawline(tree_node * &T, int flag)
{
	Scalar gray = Scalar(200, 200, 200);
	Scalar red = Scalar(100, 100, 255);
	Mat picture(500, 500, CV_8UC3, Scalar(255, 255, 255, 0.5));
	for (int i = 0; i < test_xy.size(); i++)
	{
		Point center = Point(test_xy[i][0] * 150 + 250, test_xy[i][1] * 150 + 250);
		circle(picture, center, 3, Scalar(0, 0, 0), -1);
		Point center2 = Point(test_xy[i][0] * 150 + 250, compute_tree(T, test_xy[i][0],1) * 150 + 250);
		circle(picture, center2, 3, Scalar(255, 0, 0), -1);
	}
	return picture;
}

pair<tree_node*, float> genetic_programming(int iteration, int pool, int elite, float mutation_rate, int tourment_value)
{
	vector< pair<tree_node*, float>> one_generation;
	vector<tree_node*> selection_pool;
	vector< pair<tree_node*, float>> caled_selection_pool;
	vector< pair<tree_node*, float>> next_generation;
	vector<int> survived_index;
	one_generation = init_gp(pool);
	int current_iter = 0;
	while (current_iter < iteration)
	{
		selection_pool.clear(); next_generation.clear(); caled_selection_pool.clear();
		current_iter++;
		sort(one_generation.begin(), one_generation.end(), [](const pair<tree_node*, float> &left, const pair<tree_node*, float> &right)
		{
			return left.second < right.second; 
		});
		//next_generation.push_back(one_generation[0]); // The best result will 100% be preserved
		for (int i = 0; i < elite; i++) // get_other_elite, and have certain probability to mutate it
		{
			if (random_prob(mutation_rate*0.1))
			{
				mutation_number(one_generation[i].first);
			}
			next_generation.push_back(one_generation[i]);
		}

		shuffle(begin(one_generation), end(one_generation), rng); // shuffle to cross_over
		for (int i = 0; i < pool/2; i++) // get_parents_and_children
		{
			tree_node* male = one_generation[2*i].first;
			tree_node* female = one_generation[2 * i + 1].first;
			vector<tree_node*> children = cross_over_with_limit(male, female);
			selection_pool.push_back(children[0]);
			selection_pool.push_back(children[1]);
			//selection_pool.push_back(male);
			//selection_pool.push_back(female);
		}

		for (int i = 0; i < selection_pool.size(); i++) // calculate_all_fitness
		{
			if (random_prob(mutation_rate))
			{
				mutation(selection_pool[i]);
				mutation_number(selection_pool[i]);
			}
			float fitness = get_fitness_value(selection_pool[i]);
			caled_selection_pool.push_back(make_pair(selection_pool[i], fitness));
		}

		for (int i = 0; i < next_generation.size(); i++)
		{
			next_generation[i].second = get_fitness_value(next_generation[i].first);
		}

		survived_index.clear(); // tourment_selection
		while (next_generation.size()<pool)
		{
			int index = tourment_selection(tourment_value, caled_selection_pool);
			next_generation.push_back(caled_selection_pool[index]);
			survived_index.push_back(index);
		}
		one_generation.clear();
		one_generation = next_generation;
		for (int i = 0; i < caled_selection_pool.size(); i++) // delete_lost_trees
		{
			if (std::find(survived_index.begin(), survived_index.end(), i) != survived_index.end())
				continue;
			else
			{
				//delete_tree(caled_selection_pool[i].first);
			}
		}
		if (current_iter % record_frequncy == 0) // see_result
		{
			sort(one_generation.begin(), one_generation.end(), [](const pair<tree_node*, float> &left, const pair<tree_node*, float> &right)
			{
				return left.second < right.second;
			});
			one_generation[0].second = get_fitness_value(one_generation[0].first, 0);
			if (one_generation[0].second <= best.second)
			{
				delete_tree(best.first);
				tree_node* new_best = new tree_node();
				new_best = clone_a_tree(one_generation[0].first);
				best.first = new_best;
				best.second = one_generation[0].second;
			}
			performance.push_back(make_pair(current_iter*pool, best.second));
			if (current_iter % 1000 == 0)
			{
				cout << "number_iter" << current_iter << endl;
				print_tree(best.first, "print");
				cout << endl;
				cout << "best_fitness_now: " << get_fitness_value(best.first, 0) << endl;
				Mat picture = drawline(best.first, 0);
				cv::namedWindow("demonstration");
				cv::imshow("demonstration", picture);
				cv::waitKey(1);
			}
		}
	}
	sort(one_generation.begin(), one_generation.end(), [](const pair<tree_node*, float> &left, const pair<tree_node*, float> &right)
	{
		return left.second < right.second;
	});
	if (one_generation[0].second < best.second)
	{
		best = one_generation[0];
	}
	return best;
}

void complx_and_acurcy(tree_node *&T,string filename)
{
	int complex = count_a_tree(T);
	int depth = update_depth(T, 0);
	float performance = get_fitness_value(T, 0);
	ofstream file;
	file.open(filename, fstream::app);
	file << to_string(complex) << "-" << to_string(depth) << "-" << to_string(performance) << "\n";
	file.close();
}

pair<tree_node*, float> hill_climber(int iteration, int pool)
{
	pair<tree_node*, float> best;
	vector< pair<tree_node*, float>> caled_neighbours;
	best = init_gp(1)[0];

	int current_iter = 0;
	while (current_iter < iteration)
	{
		current_iter++;
		caled_neighbours.push_back(best);
		for (int i = 0; i < pool - 1; i++)
		{
			tree_node* neighbour = new tree_node();
			clone_a_tree(neighbour);
			if (random_boolean())
				mutation(neighbour);
			else
				mutation_number(neighbour);
			float fitness = get_fitness_value(neighbour);
			caled_neighbours.push_back(make_pair(neighbour, fitness));
		}
		sort(caled_neighbours.begin(), caled_neighbours.end(), [](const pair<tree_node*, float> &left, const pair<tree_node*, float> &right)
		{
			return left.second < right.second;
		});
		best = caled_neighbours[0];
		for (int i = 1; i < pool - 1; i++)
		{
			delete_tree(caled_neighbours[i].first);
		}
		caled_neighbours.clear();
		if (current_iter % record_frequncy == 0)
		{
			performance.push_back(make_pair(current_iter*pool, best.second));
			if (current_iter % 1000 == 0)
			{
				cout << "number_iter" << current_iter << endl;
				print_tree(best.first, "print");
				cout << endl;
				cout << "best_fitness_now: " << get_fitness_value(best.first, 0) << endl;
				Mat picture = drawline(best.first, 0);
				cv::namedWindow("demonstration");
				cv::imshow("demonstration", picture);
				cv::waitKey(1);
			}
		}
	}
	return best;
}

pair<tree_node*, float> random_search(int iteration, int batch)
{
	pair<tree_node*, float> best;
	pair<tree_node*, float> candidate;
	best = init_gp(1)[0];
	int current_iter = 0;
	while (current_iter < iteration)
	{
		current_iter++;
		for (int i = 0; i < batch; i++)
		{
			candidate = init_gp(1)[0];
			if (best.second <= candidate.second)
				delete_tree(candidate.first);
			else
			{
				delete_tree(best.first);
				best = candidate;
			}
		}
		if ((current_iter) % record_frequncy == 0)
		{
			performance.push_back(make_pair(current_iter*batch, best.second));
			if ((current_iter) % 3000 == 0)
			{
				cout << "number_iter" << current_iter*30 << endl;
				print_tree(best.first, "print");
				cout << endl;
				cout << "best_fitness_now: " << get_fitness_value(best.first, 0) << endl;
				Mat picture = drawline(best.first, 0);
				cv::namedWindow("demonstration");
				cv::imshow("demonstration", picture);
				cv::waitKey(1);
			}
		}
	}
	return best;
}

void calculate_mean(string filename)
{
	float mean = 0;
	float fitness = 0;
	for (int i = 0; i < test_xy.size(); i++)
	{
		mean += test_xy[i][1];
	}
	mean = mean / test_xy.size();
	for (int i = 0; i < test_xy.size(); i++)
	{
		fitness += abs(test_xy[i][1] - mean);
	}
	ofstream file;
	file.open(filename, fstream::app);
	file << to_string(1) << "-" << to_string(1) << "-" << to_string(fitness) << "\n";
	file.close();
}

int main()
{	
	string location = "C://Users//xzy19//Desktop//GP//";
	best.second = 100;
	read_in_files("SR_div_noise.txt");
	int total_test_number = test_xy.size();
	cout << "total_test_number:" << total_test_number << endl;
	cout << "initiate_complete!" << endl << endl;
	//maxim_iteration-number -- each_generation---elite_number--mutation_rate--tourmant_number
	double t = cv::getTickCount();
	if (maxim_layer > 1)
	{
		//best = genetic_programming(10000, 30, 2, 0.1, 2);
		//best = hill_climber(10000, 30);
		best = random_search(10000,30);
		t = ((double)getTickCount() - t) / getTickFrequency();
		write_performance(location + "result_store.txt", best.second, t);
		complx_and_acurcy(best.first, location + "complex_accuracy.txt");
	}
	else
	{
		calculate_mean(location + "complex_accuracy.txt");
	}
	//cv::waitKey(1000);
	//cv::waitKey(-1);
}
