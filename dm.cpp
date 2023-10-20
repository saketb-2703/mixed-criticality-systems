#include <bits/stdc++.h>
using namespace std;

#define MAX_TASK 10

int num_of_task;
int hyperperiod;
int job_index = 0;
int preemptions = 0;

int arrival_time[MAX_TASK], period[MAX_TASK], execution_time[MAX_TASK], deadline[MAX_TASK];
vector<int> cur_task_info(MAX_TASK, 0);

struct job{
	int arrival_time;
	int period;
	int execution_time;
	int deadline;
	int task_id = -1;
	int job_id = -1;
	int job_index;
	int remain_time = 0;
	int relative_deadline;
	int cur_time;
	int preempted_time;

	bool operator<(const job& other) const {
		if(relative_deadline == other.relative_deadline)
			return preempted_time < other.preempted_time;  // The latest preempted/cpu-executed job should be selected among equal priority jobs
        return relative_deadline > other.relative_deadline; 
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
	file = fopen("tasks.txt", "r");
	if(file == NULL){
		cout << "File not found" << endl;
		exit(0);
	}
	fscanf(file, "%d", &num_of_task);
	for(int i = 0; i < num_of_task; i++){
		fscanf(file, "%d %d %d %d", &arrival_time[i], &period[i], &execution_time[i], &deadline[i]);	
	}
	hyperperiod = findLCM();
}

JOB loadJob(int task_id, int cur_time){
	JOB t;
	t.arrival_time = cur_time;
	t.period = period[task_id];
	t.execution_time = execution_time[task_id];
	t.deadline = cur_time + deadline[task_id];
	t.task_id = task_id;
	t.job_id = cur_task_info[task_id];
	t.job_index = job_index++;
	t.remain_time = execution_time[task_id];
	t.relative_deadline = deadline[task_id];
	t.cur_time = cur_time;
	t.preempted_time = -1;
	return t;
}
int findMinTask(){
	int minTime = INT_MAX;
	for(int i = 0; i < num_of_task; i++){
		minTime = min(minTime, arrival_time[i] + cur_task_info[i] * period[i]);
	}
	return minTime;
}
void deadline_monotonic(int time) {
    float utilization = 0;
	for (int i = 0; i < num_of_task; i++){
		utilization += (1.0*execution_time[i])/deadline[i];
	}
	int n = num_of_task;
	// if (utilization > 1){
	// 	cout << endl << "Given problem is not schedulable under said scheduling algorithm." << endl;
	// 	// exit(0);
	// }


    priority_queue<JOB> pq; // Min heap of pairs ===> {relative_deadline, {remain_time, task_id}}

	JOB cur_job, prev_job;
	int cur_time = 0, prev_time = 0;
    int cur_job_index = -1, prev_job_index = -1; // no job in CPU

	int decision_point = findMinTask();
    if(decision_point != cur_time)
	    printf("%d -> %d => IDLE\n", cur_time, decision_point);
	cur_time = decision_point;
	bool arrival = false; // if cur_job is going to miss deadline, we restart iteration at cur_time again. In order not to add the arriving periodic job(at this cur_time) twice, we use this flag.

    while(cur_time < time){
		if(cur_job_index != -1){
            cur_job.remain_time -= cur_time - prev_time;
            if(cur_job.remain_time == 0){  // remove completed job from CPU
                cur_job_index = -1;
				prev_job_index = -1;
			}
        }
		if(!arrival){
			for (int j = 0; j < num_of_task; j++) {
				if ((cur_time - arrival_time[j]) % period[j] == 0) {
					JOB t = loadJob(j, cur_time);
					pq.push(t);
					cur_task_info[j]++;
					if(cur_job_index != -1){
						cur_job.preempted_time = cur_time;
						pq.push(cur_job);
						cur_job_index = -1;
					}
				}
			}
			arrival = true;
		}
		if(!pq.empty()){
			cur_job = pq.top();
			pq.pop();
			cur_job_index = cur_job.job_index;

			if(cur_job.deadline <= cur_time){
                printf("Job %d of Task %d missed deadline\n", cur_job.job_id + 1, cur_job.task_id + 1);
                cur_job_index = -1;
                continue;
            }
			if(cur_job.deadline - cur_time < cur_job.remain_time){
                printf("Job %d of Task %d will miss deadline\n", cur_job.job_id + 1, cur_job.task_id + 1);
                cur_job_index = -1;
                continue;
            }

			if(prev_job_index != -1 && cur_job_index != prev_job_index){
				preemptions++;
                printf("Preemption at t = %d\n", cur_time);

			}
			prev_job_index = cur_job_index;
			
			decision_point = min(cur_time + cur_job.remain_time, findMinTask());
			
			if(decision_point == cur_time + cur_job.remain_time){
				printf("%d -> %d => T%d(J%d)*\n", cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
			}
			else{
				printf("%d -> %d => T%d(J%d)\n", cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
			}
		}
		else{
			decision_point = findMinTask();
			printf("%d -> %d => IDLE\n", cur_time, decision_point);
		}

		prev_time = cur_time;
		cur_time = decision_point;
		arrival = false;
    }
}

int main(int argc, char* argv[]){
    get_task_info();		//collecting tasks detail
	printf("Hyperperiod: %d\n", hyperperiod);
    deadline_monotonic(hyperperiod);
	printf("Preemptions: %d\n", preemptions);
	return 0;
}




// 1. Include Hyperperiod ----- DONE
// 2. Calculate missing deadline jobs prior to execution and not execute them ----- DONE
// 3. Latest/currently executing job should have more priority among equal-priority jobs ----- DONE
// 4. Include Arrival times of tasks --- DONE
// 5. Take File input ----- DONE
// 6. Struct for Jobs ----- DONE
// 7. Don't count preemption if same job selected as new job ----- DONE
// 8. Re-initialize cpu variables at actual time only and not during decision point calculation ----- DONE