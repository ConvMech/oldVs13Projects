//
//  main.cpp
//  KnapsackProblem
//
//  Created by 曹袭亚 on 2017/11/12.
//  Copyright © 2017年 曹袭亚. All rights reserved.
//
#include <vector>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;
random_device rd;
mt19937 rng(rd());

class generate_dataset{
    int object_num = 0;
    int W = 0;
    int maxim_range = 0;
    int min_range = 0;
    int min_W;
    int max_W;
    bool auto_W = true;
    vector<int> value;
    vector<int> weight;
    int random_value(int min, int max) {
        uniform_int_distribution<int>uni(min, max - 1);
        int value = uni(rng);
        return value;
    }
public:
    // number+max_v+min_v+min_W+max_W+autobool+Wifnot
    generate_dataset(int in_object_num,
                     int in_min_range, int in_max_range, int in_min_W, int in_max_W,
                     bool in_auto_W, int in_W = 0){
        object_num = in_object_num;
        maxim_range = in_max_range;
        auto_W = in_auto_W;
        min_range = in_min_range;
        W = in_W;
        min_W = in_min_W;
        max_W = in_max_W;
    }
    void random_generate(){
        for (int i = 0; i < object_num; i++){
            value.push_back(random_value(min_range, maxim_range));
            weight.push_back(random_value(min_W, max_W));
        }
        if (auto_W){
            int sum = 0;
            for (int i = 0; i < weight.size(); i++){
                sum += weight[i];
            }
            sum *= 0.75;
            W = sum;
        }
    }
    void write_data_to_file(string name){
        ofstream c;
        string cfile = name + "_c.txt";
        c.open(cfile.c_str());
        if (c.is_open()){
            c << W << endl;
            c.close();
        }
        ofstream p;
        string pfile = name + "_p.txt";
        p.open(pfile.c_str());
        if (p.is_open()){
            for (int i = 0; i < value.size(); i++){
                p << value[i] << endl;
            }
            p.close();
        }
        ofstream w;
        string wfile = name + "_w.txt";
        w.open(wfile.c_str());
        if (w.is_open()){
            for (int i = 0; i < weight.size(); i++){
                w << weight[i] << endl;
            }
            w.close();
        }
    }
};




class CItem {
    double value;
    double weight;
public:
    CItem (double v, double w) {
        value = v;
        weight = w;
    }
    double Weight(){
        return weight;
    };
    double Value(){
        return value;
    }
};

class CSolution {
    vector <bool> indicator;
    int num;
    double fitness = 0;
public:
    CSolution (int n, bool initial = true) {
        num = n;
        if (initial) {
            for (int i = 0; i < n; i++) {
                bool rdn = randomIndicator();
                indicator.push_back(rdn);
            }
        }
    }
    int SolutionSize() {
        return num;
    }
    bool randomIndicator() {
        uniform_int_distribution<int>uni(0,199);
        bool ind = (int)(uni(rng)/100);
        return ind;
    }
    double Fitness(vector<CItem> items) {
        double fit = 0;
        
        for (int i = 0; i < num; i++) {
            fit = fit + indicator[i] * items[i].Value();
        }
        fitness = fit;
        return fit;
    }
    double getFit(){
        return fitness;
    }
    void print_indicator() {
        for (int i = 0; i < num; i++)
            cout << indicator[i] << " ";
        cout << endl;
    }
    bool Constrain(vector<CItem> items, double W) {
        double c = 0;
        for (int i = 0; i < num; i++)
            c = c + indicator[i] * items[i].Weight();
        if (c < W)
            return true;
        else
            return false;
    }
    void SwapIndicator(CSolution * s1, CSolution * s2, int position) {
        int n = s1->SolutionSize();
        for (int i = 0; i < position; i++)
            indicator.push_back(s1->indicator[i]);
        for (int i = position; i < n; i++)
            indicator.push_back(s2->indicator[i]);
    }
    
    void MutateIndicator(CSolution * s, int position) {
        int n = s->SolutionSize();
        //cout << "position" << position << endl;
        int i;
        for (i = 0; i < position; i++)
            indicator.push_back(s->indicator[i]);
        //cout << "position i" << i << endl;
        if (s->indicator[position])
            indicator.push_back(false);
        else
            indicator.push_back(true);
        for (int i = position + 1; i < n; i++)
            indicator.push_back(s->indicator[i]);
    }
    
    int randomPosition(int n) {
        uniform_int_distribution<int>uni(0,n-1);
        int pos = uni(rng);
        return pos;
    }
};

CSolution * valid_init(int num, vector<CItem> items, double W)
{
	bool valid_init = false;
	CSolution * true_init = new CSolution(num);
	CSolution * init = new CSolution(num);
	valid_init = init->Constrain(items, W);
	if (valid_init)
	{
		true_init = init;
		//delete init;
		return true_init;
	}
	while (!valid_init)
	{
		//delete init;
		CSolution * init = new CSolution(num);
		valid_init = init->Constrain(items, W);
		if (valid_init)
			true_init = init;
	}
	return true_init;
}

class CPopulation {
    vector<CSolution *> population;
    int size;
    int num;
public:
    CPopulation (int s, int n, vector<CItem> items, double W, bool initial = true) {
        size = s;
        num = n;
        if (initial) {
            for (int i = 0; i < s; i++){
                bool constrain = false;
                while (!constrain) {
                    CSolution * solution = new CSolution(num);
                    constrain = solution->Constrain( items, W);
                    if (constrain){
                        solution->Fitness(items);
                        population.push_back(solution);
                        //cout << solution->SolutionSize() << endl;
                    }
                    else
                        delete solution;
                }
            }
        }
    }
    int pop_size() {
        return (int)population.size();
    }
    
    CSolution * Best(vector<CItem> items, int start, int end, bool all) {
        int best = 0;
        int best_fit = 0;
        if (end > start) {
            for (int i = start; i < end; i++) {
                int temp = population[i]->Fitness(items);
                if (temp > best_fit) {
                    best_fit = temp;
                    best = i;
                }
            }
        }
        return population[best];
    }
    
    CSolution * Specific(int position) {
        return population[position];
    }
    
    void PutIn(CSolution * s) {
        population.push_back(s);
        size++;
    }
    
    CSolution * PutOut() {
        CSolution * s = population[size - 1];
        population.pop_back();
        size--;
        return s;
    }
    
    void Sort() {
        sort(population.begin( ), population.end( ), [ ](  CSolution * s1,  CSolution * s2 ){
            return  s1->getFit() < s2->getFit();
        });
    }
    
	void shuffle_population(){
		shuffle(begin(population),end(population),rng);
	}
};

class CGeneticAlgorithm {
    double MutationRate;
    double SelectPercent;
    double CrossOverPercent;
    double ElitePercent;
public:
    CGeneticAlgorithm(double m,double s, double c, double e) {
        MutationRate = m;
        SelectPercent = s;
        CrossOverPercent = c;
        ElitePercent = e;
    }
    vector<CSolution *> CrossOver(CSolution * s1, CSolution * s2) {
        int n = s1->SolutionSize();
        CSolution * new_s1 = new CSolution(n, false);
        CSolution * new_s2 = new CSolution(n, false);
        int pos = s1->randomPosition(n);
        new_s1->SwapIndicator(s1, s2, pos);
        new_s2->SwapIndicator(s2, s1, pos);
        vector<CSolution *> new_s;
        new_s.push_back(new_s1);
        new_s.push_back(new_s2);
        return new_s;
    }
	vector<CSolution *> Valid_CrossOver(CSolution * s1, CSolution * s2, vector<CItem> items, double W) {
		int childern_number = 2;
		vector<CSolution *> new_s;
		int i = 0;
		while (i < childern_number)
		{
			vector<CSolution *> childs = CrossOver(s1, s2);
			if (childs[0]->Constrain(items, W))
			{
				i++;
				new_s.push_back(childs[0]);
			}
			if (childs[1]->Constrain(items, W))
			{
				i++;
				new_s.push_back(childs[1]);
			}
		}
		return new_s;
	}
    
    CSolution * SelectParent(CPopulation POP, vector<CItem> items) {
        int pop_size = POP.pop_size();
        CSolution temp(0, false);
        int s = temp.randomPosition(pop_size);
        int e = int(s + pop_size * SelectPercent) % pop_size;
        CSolution * parent = POP.Best(items, s, e, false);
        return parent;
    }
    
    CSolution * Mutation(CSolution * s) {
        int n = s->SolutionSize();
        int mutatePos = s->randomPosition(n);
        CSolution * new_s = new CSolution(n,false);
        new_s->MutateIndicator(s, mutatePos);
        return new_s;
    }

	CPopulation Evolve_T(CPopulation POP, vector<CItem> items, double W){
		int size = POP.pop_size();
		int n = (int)items.size();
		int co_size = size * CrossOverPercent;
		int mu_size = size * MutationRate;
		int el_size = size * ElitePercent;
		CPopulation selection_pool(0, n, items, W, false);
		CPopulation mutated_pool(0, n, items, W, false);
		CPopulation elite_pop(0, n, items, W, false);
		CPopulation next_gen(0, n, items, W, false);

		POP.Sort(); // ready_to_select_elite;
		for (int i = 0; i < el_size; i++) {
			CSolution * s = POP.PutOut();
			elite_pop.PutIn(s);
			//cout << "elite: " << s->getFit() << endl;
		}
		POP.shuffle_population(); // shuffle the population to cross_over;
		while (POP.pop_size() >0){
			CSolution * s1 = POP.PutOut();
			CSolution * s2 = POP.PutOut();
			vector<CSolution *> children = Valid_CrossOver(s1, s2, items, W);
			children[0]->Fitness(items);
			selection_pool.PutIn(children[0]);
			children[1]->Fitness(items);
			selection_pool.PutIn(children[1]);
		}
		selection_pool.shuffle_population();

		vector <int >two_container;
		for (int i = 0; i < selection_pool.pop_size(); i ++)
		{
			two_container.push_back(2 * i);
			two_container.push_back(2 * i + 1);
			CSolution * selected = selection_pool.Specific(i);
			CSolution * reborn = Mutation(selected);
			if (reborn->Constrain(items, W)){
				reborn->Fitness(items);
				mutated_pool.PutIn(reborn);
			}
			else{
				CSolution * reborn_ = valid_init(items.size(), items, W);
				reborn_->Fitness(items);
				mutated_pool.PutIn(reborn_);
			}
		}

		shuffle(begin(two_container), end(two_container), rng);
		
		for (int i = 0; i < selection_pool.pop_size(); i++)
		{
			CSolution *player1 = selection_pool.Specific(two_container[2 * i]%selection_pool.pop_size());
			CSolution *player2 = mutated_pool.Specific(two_container[2 * i + 1] % selection_pool.pop_size());
			if (player1->getFit() > player2->getFit())
				next_gen.PutIn(player1);
			else
				next_gen.PutIn(player2);
		}

		for (int i = 0; i < elite_pop.pop_size(); i++) {
			next_gen.PutIn(elite_pop.Specific(i));
		}

	/*	cout << size << endl;
		cout << elite_pop.pop_size() << endl;
		cout << selection_pool.pop_size()<< endl;
		cout << mutated_pool.pop_size() << endl;
		cout << next_gen.pop_size() << endl;*/

		return next_gen;
	}
    
    CPopulation Evolve(CPopulation POP, vector<CItem> items, double W) {
        int size = POP.pop_size();
        int n = (int)items.size();
        int co_size = size * CrossOverPercent;
        int mu_size = size * MutationRate;
        int el_size = size * ElitePercent;
        CPopulation crossover_pop(0, n, items, W, false);
        CPopulation mutate_pop(0, n, items, W, false);
        CPopulation elite_pop(0, n, items, W, false);
        
        
        for (int i = 0; i < co_size; ) {
            CSolution * s1 = SelectParent(POP, items);
            CSolution * s2 = SelectParent(POP, items);
            vector<CSolution *> childs = CrossOver(s1, s2);
            if (childs[0]->Constrain(items, W)) {
                childs[0]->Fitness(items);
                //cout << "child "<<childs[0]->getFit()<<endl;
                crossover_pop.PutIn(childs[0]);
                i++;
            }
            if (childs[1]->Constrain(items, W)) {
                childs[1]->Fitness(items);
                crossover_pop.PutIn(childs[1]);
                i++;
            }
        }
        
        for (int i = 0; i < mu_size; ) {
            CSolution temp(0, false);
            int pos = temp.randomPosition(size);
            CSolution * s = POP.Specific(pos);
            //cout << "pos: " << pos << endl;
            //cout << "solution size: " << s->SolutionSize() << endl;
            CSolution * new_s = Mutation(s);
            if (new_s->Constrain(items, W)){
                i++;
                new_s->Fitness(items);
                mutate_pop.PutIn(new_s);
            }
        }
        
        POP.Sort();
        for (int i = 0; i < el_size; i++) {
            CSolution * s = POP.PutOut();
            elite_pop.PutIn(s);
            //cout << "elite: " << s->getFit() << endl;
        }
                
        for (int i = el_size; i < size; i++) {
            CSolution * s = POP.PutOut();
            delete s;
        }
        //cout << "co_pop size: " << crossover_pop.pop_size() << endl;
        //cout << "mu_pop size: " << mutate_pop.pop_size() << endl;
        CPopulation new_pop = SelectionMethod(size, elite_pop, crossover_pop, mutate_pop, items, W);
        return new_pop;
    }
    
    CPopulation SelectionMethod(int size, CPopulation elite_pop, CPopulation co_pop, CPopulation mu_pop, vector<CItem> items, double W) {
        int n = (int)items.size();
        CPopulation new_pop(0, n, items, W, false);
        for (int i = 0; i < elite_pop.pop_size(); i++) {
            new_pop.PutIn(elite_pop.Specific(i));
        }
        co_pop.Sort();
        mu_pop.Sort();
        for (int i = elite_pop.pop_size(); i < size; i++) {
            //cout << "co_pop size: " << co_pop.pop_size() << endl;
            //cout << "mu_pop size: " << mu_pop.pop_size() << endl;
            CSolution * elite1 = co_pop.PutOut();
            CSolution * elite2 = mu_pop.PutOut();
            if (elite1->getFit() > elite2->getFit())
                new_pop.PutIn(elite1);
            else
                new_pop.PutIn(elite1);
        }
        //cout << "new_pop: " << new_pop.pop_size() << endl;
        return new_pop;
    }
};

class CHill_climber
{
    int iteration_limit = 0;
    int num;
    vector<CItem> items;
    double W;
public:
    CHill_climber(int iteration_limit_input, int num_input, vector<CItem> items_input, double W_input)
    {
        iteration_limit = iteration_limit_input / num_input;
        cout << iteration_limit_input << " " <<num_input << "  test"<<endl;
        num = num_input;
        items = items_input;
        W = W_input;
    }
    CSolution * climb_one_step(CSolution* origin)
    {
        CSolution *seed = origin;
        CPopulation step_neibours(0, num, items, W, false);
        step_neibours.PutIn(seed);
        for (int i = 0; i < num; i++)
        {
            CSolution *one_neighbour = new CSolution(num, false);
            one_neighbour->MutateIndicator(seed, i);
            
            if (one_neighbour->Constrain(items, W))
            {
                one_neighbour->Fitness(items);
                step_neibours.PutIn(one_neighbour);
            }
            else
            {
                bool valid_sub= false;
                CSolution * true_sub = new CSolution(num);
                CSolution * substitue= new CSolution(num);
                valid_sub = substitue->Constrain(items, W);
                if (valid_sub)
                    true_sub = substitue;
                else
                    delete substitue;
                while (!valid_sub)
                {
                    CSolution * substitue = new CSolution(num);
                    valid_sub = substitue->Constrain(items, W);
                    if (valid_sub)
                        true_sub = substitue;
                    else
                        delete substitue;
                }
                true_sub->Fitness(items);
                step_neibours.PutIn(true_sub);
            }
        }
        //cout << step_neibours.pop_size() << endl;
        step_neibours.Sort();
        CSolution * new_origin = step_neibours.PutOut();
        while(step_neibours.pop_size() > 0){
            CSolution * new_origin = step_neibours.PutOut();
            delete new_origin;
        }
        //cout << new_origin->getFit() << endl;
        return new_origin;
    }
    CSolution * climb(string problem, int per_rec)
    {
        string resultfile = problem + "_hc_result.txt";
        ofstream outfile(resultfile, ios::app);
        outfile << resultfile << endl;
        bool valid_init = false;
        CSolution * true_init = new CSolution(num);
        CSolution * init = new CSolution(num);
        valid_init = init->Constrain(items, W);
        if (valid_init)
            true_init = init;
        else
        {
            while (!valid_init){
                delete init;
                CSolution * init = new CSolution(num);
                valid_init = init->Constrain(items, W);
                if (valid_init)
                    true_init = init;
            }
        }
        cout << (true_init->Fitness(items)) << endl;
        cout << iteration_limit << endl;
        for (int i = 0; i < iteration_limit; i++)
        {
            true_init = climb_one_step(true_init);
            //if (i % per_rec == 0)
            outfile << true_init->getFit() << endl;
            cout << true_init->getFit() << endl;
        }
        outfile.close();
        return true_init;
    }
};

class CRandom {
    int iteration;
    int num;
    vector<CItem> items;
    double W;
public:
    CRandom(vector<CItem> items_, double W_, int num_, int iteration_) {
        iteration = iteration_;
        num = num_;
        W = W_;
        items = items_;
    }
    
    CSolution * random_one_step(CSolution* origin)
    {
        bool better = false;
        int count = 0;
        while(!better){
            count++;
            CSolution * new_s= new CSolution(num, false);
            int pos = new_s->randomPosition(num);
            new_s->MutateIndicator(origin, pos);
            if (new_s->Constrain(items, W)) {
                new_s->Fitness(items);
                if (new_s->getFit() > origin->getFit()) {
                    better = true;
                    return new_s;
                }
                else
                    delete new_s;
            }
            else
                delete new_s;
            if (count > 50)
                break;
        }
        return origin;
    }
    
    CSolution * Random(string problem, int per_rec) {
        string resultfile = problem + "_rand_result.txt";
        ofstream outfile(resultfile, ios::app);
        outfile << resultfile << endl;
        bool valid_init = false;
        CSolution * true_init = new CSolution(num);
        CSolution * init = new CSolution(num);
        valid_init = init->Constrain(items, W);
        if (valid_init)
            true_init = init;
        else
            delete init;
        while (!valid_init){
            CSolution * init = new CSolution(num);
            valid_init = init->Constrain(items, W);
            if (valid_init)
                true_init = init;
            else
                delete init;
        }
        cout << (true_init->Fitness(items)) << endl;
        for (int i = 0; i < iteration; i++){
            true_init = random_one_step(true_init);
            if (i % per_rec == 0)
                outfile << true_init->getFit() << endl;
        }
        outfile.close();
        return true_init;
    }
};

class CSimulated_annealing
{
    int iteration_limit = 0;
    int num;
    vector<CItem> items;
    double W;
public:
    CSimulated_annealing(int iteration_limit_input, int num_input, vector<CItem> items_input, double W_input)
    {
        iteration_limit = iteration_limit_input;
        num = num_input;
        items = items_input;
        W = W_input;
    }
    
    CSolution * init()
    {
        bool valid_init = false;
        CSolution * true_init = new CSolution(num);
        CSolution * init = new CSolution(num);
        valid_init = init->Constrain(items, W);
        if (valid_init)
            return true_init = init;
        while (!valid_init)
        {
            delete init;
            CSolution * init = new CSolution(num);
            valid_init = init->Constrain(items, W);
            if (valid_init)
                true_init = init;
        }
        return true_init;
    }
    CSolution * mutate(CSolution * mutate_base)
    {
        int n = mutate_base->SolutionSize();
        CSolution temp(0, false);
        int mutatePos = temp.randomPosition(n);
        CSolution * new_s = new CSolution(n, false);
        new_s->MutateIndicator(mutate_base, mutatePos);
        
        if (new_s->Constrain(items, W))
            return new_s;
        else
            delete new_s;
        bool valid = false;
        int count = 0;
        while (!valid)
        {
            count++;
            if (count > 50)
                break;
            int n = mutate_base->SolutionSize();
            CSolution temp(0, false);
            int mutatePos = temp.randomPosition(n);
            CSolution * new_s = new CSolution(n, false);
            new_s->MutateIndicator(mutate_base, mutatePos);
            if (!(new_s->Constrain(items, W)))
            {
                delete new_s;
                continue;
            }
            else
                valid = true;
            return new_s;
        }
        return mutate_base;
    }
    double random_0to1()
    {
        uniform_int_distribution<int> uni(1, 100000);
        return((double)uni(rng) / 100000);
    }
    CSolution *annel(string problem, int per_rec)
    {
        string resultfile = problem + "_annel_result.txt";
        ofstream outfile(resultfile, ios::app);
        outfile << resultfile << endl;
        double completed_portion = 0;
        int iter_count = 0;
        CSolution * start= init();
        while (completed_portion < 1)
        {
            float defender = start->Fitness(items);
            iter_count++;
            completed_portion = (double)(iter_count) / iteration_limit;
            double temperatue = completed_portion;
            CSolution *mutated = mutate(start);
            float candidate = mutated->Fitness(items);
            if (candidate >= defender)
                start = mutated;
            else
            {
                float probability = 0.6 * ((candidate - defender) / defender) + 1.3*temperatue;
                if (probability < random_0to1())
                    start = mutated;
                else
                {
                    delete mutated;
                    start = start;
                }
            }
            if (iter_count % per_rec == 0)
                outfile << start->getFit() << endl;
            //cout<<start->getFit()<<endl;
        }
        return start;
    }
};
class CRun {
    int num;
    vector<CItem> items;
    double W;
    string problem;
public:
    CRun(string prefix) {
        problem = prefix;
        items = ReadData(prefix);
        W = ReadConstrain(prefix);
        num = (int)items.size();
    }
    void ga(int popsize, int round, double mu_rate, double select_p, double co_rate, double eli_rate) {
        string resultfile = problem + "_ga_result.txt";
        ofstream outfile(resultfile, ios::app);
        outfile << resultfile << endl;
        CPopulation generation(popsize, num, items, W);
        CGeneticAlgorithm GA(mu_rate, select_p, co_rate, eli_rate);
        for (int i = 0; i < round; i++) {
            //generation = GA.Evolve(generation, items, W);
			generation = GA.Evolve_T(generation, items, W);
            generation.Sort();
            CSolution * temp = generation.PutOut();
			if ( i%10 == 0)
				cout << "Best" << temp->getFit() << endl;
            outfile << temp->getFit() << endl;
            generation.PutIn(temp);
        }
        outfile.close();
        generation.Sort();
        (generation.PutOut())->print_indicator();
        cout << "Hello, World!\n";
    }
    
    void hc(int iteration, int per_rec) {
        CHill_climber hill_climber(iteration, num, items, W);
        CSolution * result = hill_climber.climb(problem, per_rec);
        result->print_indicator();
        cout << result->Fitness(items) << endl;
    }
    
    void random(int iteration, int per_rec) {
        CRandom rand(items, W, num, iteration);
        CSolution * result = rand.Random(problem, per_rec);
        result->print_indicator();
        cout << result->Fitness(items) << endl;
        cout << "Hello, World!\n";
    }
    
    void sa(int iteration, int per_rec) {
        CSimulated_annealing simulated_annealing(iteration, num, items, W);
        CSolution * result = simulated_annealing.annel(problem, per_rec);
        result->print_indicator();
        cout << result->Fitness(items) << endl;
        cout << "Hello, World!\n";
    }
    double ReadConstrain(string prefix) {
        double W;
        string filename;
        filename = prefix + "_c.txt";
        fstream in(filename);
        in >> W;
        //cout << "Weight: "<< W << endl;
        return W;
    }
    vector <CItem> ReadData(string prefix) {
        vector<CItem> items;
        string vfilename, wfilename;
        vfilename = prefix + "_p.txt";
        wfilename = prefix + "_w.txt";
        fstream vin(vfilename);
        fstream win(wfilename);
        double value, weight;
        while (!vin.eof()){
            vin >> value;
            win >> weight;
            CItem item(value, weight);
            items.push_back(item);
            //cout << "value: " << value << " weight: " << weight << endl;
        }
        return items;
    }
};

int main(int argc, const char * argv[]) {
    
    string prefix = "/Users/caoxiya/Desktop/D/project/KnapsackProblem/dataset/test";
    CRun run(prefix);
    int popsize = 100;
    int round = 1e3;
    int eval = popsize * round;
    cout << "thanksgiving" << endl;
    for (int i = 0; i < 1; i++) {
        clock_t start, end;
        start = clock();
        run.ga(popsize, round, 1.2, 0.3, 1.2, 0.02);
        end = clock();
        cout << "Total Time: " <<  (double)(end - start) / CLOCKS_PER_SEC<< endl;
        start = clock();
        run.hc(eval, popsize);
        end = clock();
        cout << "Total Time: " <<  (double)(end - start) / CLOCKS_PER_SEC<< endl;
        start = clock();
        run.random(eval, popsize);
        end = clock();
        cout << "Total Time: " <<  (double)(end - start) / CLOCKS_PER_SEC<< endl;
        start = clock();
        run.sa(eval, popsize);
        end = clock();
        cout << "Total Time: " <<  (double)(end - start) / CLOCKS_PER_SEC<< endl;
    }
    cout << "thanksgiving" << endl;
    
    
    /*generate_dataset generator(1000, 1e3, 1e4, 1e2, 1e3, true);
    generator.random_generate();
    generator.write_data_to_file(prefix);*/
    return 0;
}
