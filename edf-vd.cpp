#include <bits/stdc++.h>
#include<unistd.h>  
using namespace std;

#define MAX_TASK 10
#define MAX_CRITICALITY 10

int num_of_task;
int hyperperiod;
int job_index = 0;
int preemptions = 0;
int max_criticality = 0;

double arrival_time[MAX_TASK], execution_time[MAX_TASK][MAX_CRITICALITY], deadline[MAX_TASK], virtualDeadlines[MAX_TASK];
int criticality[MAX_TASK], period[MAX_TASK];
vector<int> cur_task_info(MAX_TASK, 0);

template<typename T>
class custom_priority_queue : public priority_queue<T, vector<T>>{
    public:
        bool remove(const T& job) {
            for(auto it = this->c.begin(); it != this->c.end(); ++it) {
                if(it->job_index == job.job_index) {
                    this->c.erase(it);
                    make_heap(this->c.begin(), this->c.end(), this->comp);
                    return true;
                }
            }
            return false;
        }
};

struct job{
	double arrival_time;
	int period;
	double executed_time;
	double deadline;
    double original_deadline;
	int task_id = -1;
	int job_id = -1;
	int job_index;
	double remain_time;
	double relative_deadline;
	double cur_time;
	double preempted_time;
    int criticality;

	bool operator<(const job& other) const {
		if(deadline == other.deadline) 
			return preempted_time < other.preempted_time; // The latest preempted/cpu-executed job should be selected among equal priority jobs
        return deadline > other.deadline; 
    }
};

typedef struct job JOB;

int gcd(int a, int b) {
    if (b == 0)
        return a;
    return gcd(b, a % b);
}

int lcm(int a, int b) {
    return (a * b) / gcd(a, b);
}

// Function to calculate the LCM of an array of numbers
int findLCM() {
    int result = period[0];
    
    for (int i = 1; i < num_of_task; i++) {
        result = lcm(result, period[i]);
    }
    
    return result;
}
//collecting details of tasks
void get_task_info(){
	FILE* file;
	file = fopen("tasks-vd.txt", "r");
	if(file == NULL){
		cout << "File not found" << endl;
		exit(0);
	}
	fscanf(file, "%d", &num_of_task);
	for(int i = 0; i < num_of_task; i++){
		fscanf(file, "%lf %d %lf %d", &arrival_time[i], &period[i], &deadline[i], &criticality[i]);
        max_criticality = max(max_criticality, criticality[i]);
        for(int j = 0; j < criticality[i]; j++){
            fscanf(file, "%lf", &execution_time[i][j]);
        }
	}
    hyperperiod = findLCM();
}

JOB loadJob(int task_id, double cur_time){
	JOB t;
	t.arrival_time = cur_time;
	t.period = period[task_id];
	t.executed_time = 0;
    t.deadline = cur_time + virtualDeadlines[task_id];
    t.original_deadline = cur_time + deadline[task_id];
	t.task_id = task_id;
	t.job_id = cur_task_info[task_id];
	t.job_index = job_index++;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(0.0, execution_time[task_id][criticality[task_id]-1]);
	t.remain_time = dis(gen);

    // printf("T%d(%d) remains %lf\n", t.task_id + 1, t.job_id + 1, t.remain_time);
	t.relative_deadline = deadline[task_id];
	t.cur_time = cur_time;
	t.preempted_time = -1;
    t.criticality = criticality[task_id];
	return t;
}

double findMinTask(){
	double minTime = DBL_MAX;
	for(int i = 0; i < num_of_task; i++){
		minTime = min(minTime, arrival_time[i] + cur_task_info[i] * period[i]);
	}
	return minTime;
}
double worstUtilization(int max_criticality){
    double utilization = 0;
    // for(int i = 0; i < num_of_task; i++){
    //     utilization += (1.0*execution_time[i][criticality[i]-1])/deadline[i];
    // }
    for(int i = 1; i <= max_criticality; i++){
        for(int j = 0; j < num_of_task; j++){
            if(criticality[j] == i)
                utilization += (1.0*execution_time[j][i-1])/deadline[j];
        }
    }
    return utilization;
}
double utilisation(int l, int k){
    double utilization = 0;
    for(int i = 0; i < num_of_task; i++){
        if(criticality[i] == l){
            utilization += (1.0*execution_time[i][k-1])/deadline[i];
        }
    }
    return utilization;
}
bool checkKCondition(int k){
    double lhs_num = 0.0, lhs_den = 0.0, rhs_num = 0.0, rhs_den = 0.0;
    for(int i = k+1; i <= max_criticality; i++){
        lhs_num += utilisation(i, k);
    }
    for(int i = 1; i <= k; i++){
        lhs_den += utilisation(i, i);
    }
    lhs_den = 1 - lhs_den;
    for(int i = k+1; i <= max_criticality; i++){
        rhs_num += utilisation(i, i);
    }
    rhs_num = 1 - rhs_num;
    for(int i = 1; i <= k; i++){
        rhs_den += utilisation(i, i);
    }
    double lhs = lhs_num/lhs_den;
    double rhs = rhs_num/rhs_den;
    lhs = round(lhs * 1e6) / 1e6;
    rhs = round(rhs * 1e6) / 1e6;
    if(lhs <= rhs)
        return true;
    return false;
}
int findK(){
    vector<int> possibleK;
    for(int i = 1; i <= max_criticality; i++){
        if(1 - worstUtilization(i) <= 0)
            continue;
        if(checkKCondition(i))
            possibleK.push_back(i);
    }
    int size = possibleK.size();
    return size == 0 ? -1 : possibleK[size/2];
}
double findX(int k){
    double lhs_num = 0.0, lhs_den = 0.0;
    for(int i = k+1; i <= max_criticality; i++){
        lhs_num += utilisation(i, k);
    }
    for(int i = 1; i <= k; i++){
        lhs_den += utilisation(i, i);
    }
    lhs_den = 1 - lhs_den;
    double lhs = lhs_num/lhs_den;
    lhs = round(lhs * 1e6) / 1e6;
    return lhs;
}
int offlinePreprocessing(){
    if(worstUtilization(max_criticality) <= 1){
        for(int i = 0; i < num_of_task; i++){
            virtualDeadlines[i] = 1.0*deadline[i];
        }
        return max_criticality;
    }
    else{
        int k = findK();
        if(k == -1){
            cout << "No solution exists" << endl;
            // return 1;
            exit(0);
        }
        double x = findX(k);
        // cout << "k = " << k << " x = " << x << endl;
        for(int i = 0; i < num_of_task; i++){
            if(criticality[i] > k)
                virtualDeadlines[i] = x * deadline[i];
            else    
                virtualDeadlines[i] = 1.0*deadline[i];
        }
        return k;
    }
}
void runtimeScheduling(int k, int time){

    priority_queue<JOB> pqLO; // Min heap of pairs, where first is deadline and second is task id
    custom_priority_queue<JOB> pqHI; // Min heap of pairs, where first is deadline and second is task id

	JOB cur_job;
	double cur_time = 0.0, prev_time = 0.0;
	int cur_job_index = -1, prev_job_index = -1; // no job in CPU
    int current_level = 1;

	double decision_point = findMinTask();
    if(decision_point != cur_time)
	    printf("%lf -> %lf => IDLE\n", cur_time, decision_point);
	cur_time = decision_point;
	bool arrival = false;

    while(cur_time < time){
		if(cur_job_index != -1){
            cur_job.remain_time -= cur_time - prev_time;
            cur_job.executed_time += cur_time - prev_time;

            cur_job.executed_time = round(cur_job.executed_time * 1e6) / 1e6;
            if(cur_job.executed_time >= execution_time[cur_job.task_id][current_level-1]){ // CHANGE >= to >
                printf("\nCriticality Level changed from %d to %d\n\n", current_level, current_level+1);
                current_level++;
            }


            if(abs(cur_job.remain_time) <= 0.0000001){  // remove completed job from CPU
                // printf("Job %d of Task %d completed\n", cur_job.job_id + 1, cur_job.task_id + 1);
                cur_job_index = -1; 
				prev_job_index = -1;
			}

            
            if(cur_job_index != -1){
                // printf("Job %d of Task %d preempted(%lf)\n", cur_job.job_id + 1, cur_job.task_id + 1, cur_job.remain_time);
                cur_job.preempted_time = cur_time;  // current_level changing is also a preemption point right now
                if(cur_job.criticality <= k)
                    pqLO.push(cur_job);
                else{
                    pqLO.push(cur_job);
                    cur_job.deadline = cur_job.original_deadline;
                    pqHI.push(cur_job);
                }
                cur_job_index = -1;
            }
        }
		if(!arrival){
			for (int j = 0; j < num_of_task; j++) {
				if (fmod((cur_time - arrival_time[j]), period[j]) == 0) {
					JOB t = loadJob(j, cur_time);
					pqLO.push(t);
                    if(criticality[j] > k){
                        t.deadline = cur_time + deadline[j];
                        pqHI.push(t);
                    }
					cur_task_info[j]++;
				}
			}
			arrival = true;
		}
		if((current_level <= k && !pqLO.empty()) or (current_level > k && !pqHI.empty())){
            int flag = 0;
            if(current_level <= k){
                while(!pqLO.empty()){
                    cur_job = pqLO.top();
                    pqLO.pop();
                    if(cur_job.criticality > k){
                        bool tmp = pqHI.remove(cur_job);
                        // cout << tmp << endl;
                    }
                    if(cur_job.criticality >= current_level){
                        flag++;
                        break;
                    }
                }
            }
            else{
                while(!pqHI.empty()){
                    cur_job = pqHI.top();
                    pqHI.pop();
                    if(cur_job.criticality >= current_level){
                        flag++;
                        break; 
                    }
                }
            }
            if(!flag){ // If none of the jobs in the queue have criticality >= current_level, restart loop --- In that case IDLE will be printed in next loop's else case since queue is empty
                cur_job_index = -1;
                continue;
            }
			cur_job_index = cur_job.job_index;

			if(cur_job.original_deadline <= cur_time){
                printf("Job %d of Task %d missed deadline(%lf)\n", cur_job.job_id + 1, cur_job.task_id + 1, cur_job.deadline);
                cur_job_index = -1;
                continue;
            }
			if(cur_job.original_deadline - cur_time < cur_job.remain_time){
                printf("Job %d of Task %d will miss deadline\n", cur_job.job_id + 1, cur_job.task_id + 1);
                cur_job_index = -1;
                continue;
            }

			if(prev_job_index != -1 && cur_job_index != prev_job_index){
				preemptions++;
                printf("Preemption at t = %lf\n", cur_time);
			}
			prev_job_index = cur_job_index;

			decision_point = min({cur_time + cur_job.remain_time, cur_time + execution_time[cur_job.task_id][current_level-1] -  cur_job.executed_time, findMinTask()});

			// printf("%lf %lf %lf\n", cur_job.remain_time, execution_time[cur_job.task_id][current_level-1] -  cur_job.executed_time, findMinTask());

			if(decision_point == cur_time + cur_job.remain_time){
				printf("%lf -> %lf => T%d(J%d)*\n", cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
			}
			else{
				printf("%lf -> %lf => T%d(J%d)\n", cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
			}
		}
		else{
			decision_point = findMinTask();
			printf("%lf -> %lf => IDLE\n", cur_time, decision_point);
		}
        usleep(100000);
		prev_time = cur_time;
		cur_time = decision_point;
		arrival = false;

    }
}
void edfvd(int time) {
    int k = offlinePreprocessing();
    runtimeScheduling(k, time);
}

int main(int argc, char* argv[]){
    get_task_info();		//collecting tasks detail
	printf("Hyperperiod: %d\n", hyperperiod);
    edfvd(hyperperiod);
	printf("Preemptions: %d\n", preemptions);
	return 0;
}





/*
2
0 4 4 1 2
0 6 6 2 1 5
*/