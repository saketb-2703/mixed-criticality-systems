#include <bits/stdc++.h>
#include<unistd.h>  
using namespace std;

#define MAX_TASK 10
#define MAX_JOBS 100
#define MAX_CRITICALITY 10
#define SLEEP_TIME 10000
#define COMPLETION_LB 0.8
#define PROCRAST_THRESHOLD 0.5

int num_of_task;
int num_of_jobs;
int hyperperiod;
int job_index = 0;
int preemptions = 0;
int max_criticality = 0;
int current_level = 1;
int dynamic_discarded = 0;

double arrival_time[MAX_TASK], execution_time[MAX_TASK][MAX_CRITICALITY], deadline[MAX_TASK], virtualDeadlines[MAX_TASK], slackTable[MAX_CRITICALITY][MAX_JOBS][MAX_JOBS];
int criticality[MAX_TASK], period[MAX_TASK];
vector<int> cur_task_info(MAX_TASK, 0);

typedef struct job JOB;
vector<JOB> jobs;
vector<JOB> current_jobs(MAX_TASK);
priority_queue<JOB> discardedJobs;

double total_idle_time = 0.0, total_stolen_time = 0.0, slackSystem = 0.0;


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
	double deadline = -1;
    double original_deadline;
	int task_id = -1;
	int job_id = -1;
	int job_index;
    int deadline_index;
	double remain_time;
	double relative_deadline;
	double cur_time;
	double preempted_time = -1;
    int criticality;
    bool discarded = false;

	bool operator<(const job& other) const {
		if(deadline == other.deadline) 
			return preempted_time < other.preempted_time; // The latest preempted/cpu-executed job should be selected among equal priority jobs
        return deadline > other.deadline; 
    }
};


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

// int findPartition(JOB job, vector<double> partitionDeadlines){
//     int partition = num_of_task-1;
//     for(int i = 0; i < num_of_task-1; i++){
//         if(job.deadline >= partitionDeadlines[i] && job.deadline < partitionDeadlines[i+1]){
//             partition = i;
//             break;
//         }
//     }
//     return partition;
// }
// double findSumExeCurrentJobs(int partition, vector<JOB> current_jobs){
//     double sum = 0.0;
//     for(int i = partition; i < num_of_task; i++){
//         sum += current_jobs[i].executed_time;
//     }
//     return sum;
// }

void transferRq(priority_queue<JOB>& pqLO){
    priority_queue<JOB> tmp;
    while(!pqLO.empty()){
        JOB cur_job = pqLO.top();
        pqLO.pop();
        cur_job.deadline = cur_job.original_deadline;
        tmp.push(cur_job);
    }
    while(!tmp.empty()){
        JOB cur_job = tmp.top();
        tmp.pop();
        pqLO.push(cur_job);
    }
    // pqLO = tmp;
}

void loadJobsInHyperperiod(int hyperperiod){
    vector<int> task_info(MAX_TASK, 0);
    for (int time = 0; time < hyperperiod; time++) {
        for (int i = 0; i < MAX_TASK; i++) {
            if (fmod((time - arrival_time[i]), period[i]) == 0) {
                JOB job;
                job.task_id = i;
                job.job_id = task_info[i]++;
                job.job_index = jobs.size();
                job.arrival_time = time;
                job.deadline = time + virtualDeadlines[i];  // virtualDeadlines
                job.original_deadline = time + deadline[i];
                job.remain_time = execution_time[i][1]; // at base criticality
                job.criticality = criticality[i];
                jobs.push_back(job);
            }
        }
    }
    sort(jobs.begin(), jobs.end(), [](JOB a, JOB b) {
        return a.deadline < b.deadline;
    });
    for(int i = 0; i < jobs.size(); i++){
        jobs[i].deadline_index = i;
    }
}


void slackTablePreCompute(){
    vector<vector<double>> initialSlack(max_criticality, vector<double>(num_of_jobs));
    for(int p = 1; p <= max_criticality; p++){
        for (int i = 0; i < num_of_jobs; i++) {
            initialSlack[p-1][i] = jobs[i].deadline - 0; // Assuming time origin is 0
            for (int j = 0; j < num_of_jobs; j++) {
                if((jobs[j].criticality >= p) && (jobs[j].deadline <= jobs[i].deadline))  // need to consider execution times of jobs with criticality >= p[because other jobs are discarded] and deadline less than or equal to deadline of job i
                    initialSlack[p-1][i] -= execution_time[jobs[j].task_id][p-1];
            }
            // cout << initialSlack[p-1][i] << ",";
        }
        // cout << endl;
    }
    for(int p = 1; p <= max_criticality; p++){
        for(int i = 0; i < num_of_jobs; i++){
            double miniSlack = INT_MAX;
            for(int j = i; j < num_of_jobs; j++){
                miniSlack = min(miniSlack, initialSlack[p-1][j]);
                slackTable[p-1][i][j] = miniSlack;
            }
        }
    }
}

int handleArrivals(double cur_time, int k, priority_queue<JOB>& pqLO){
    int discardedCount = 0;
    for (int j = 0; j < num_of_task; j++) {
        if (fmod((cur_time - arrival_time[j]), period[j]) == 0) {
            JOB t = loadJob(j, cur_time);
            if(t.criticality >= current_level)
                pqLO.push(t);
            else{
                t.discarded = true;
                discardedJobs.push(t);
                discardedCount++;
            }
            // printf("Job %d of Task %d arrived with criticality %d\n", t.job_id + 1, t.task_id + 1, t.criticality);
            cur_task_info[j]++;
        }
    }
    return discardedCount;
      
}

// double findPreSlack(int left, int right, double& deadline){
//     vector<vector<double>> initialSlack(max_criticality, vector<double>(num_of_jobs));
//     double miniSlack = INT_MAX;

//     for (int i = left; i <= right; i++) {
//         initialSlack[current_level-1][i] = min(deadline, jobs[i].deadline) - 0; // Assuming time origin is 0
//         for (int j = 0; j < num_of_jobs; j++) {
//             if(((jobs[j].criticality >= current_level) or (i == j)) && (jobs[j].deadline <= deadline))  // need to consider execution times of jobs with criticality >= p[because other jobs are discarded] and deadline less than or equal to deadline of job i
//                 initialSlack[current_level-1][i] -= execution_time[jobs[j].task_id][current_level-1];
//             else if(((jobs[j].criticality >= current_level) or (i == j)) && (jobs[j].deadline <= jobs[i].deadline)){ 
//                 double exeFrac = (deadline - jobs[j].arrival_time)/(jobs[j].deadline - jobs[j].arrival_time);
//                 initialSlack[current_level-1][i] -= exeFrac * execution_time[jobs[j].task_id][current_level-1];
//             }
//         }
//         miniSlack = min(miniSlack, initialSlack[current_level-1][i]);
//     }

//     return miniSlack;
// }
// void computeCurrentSlack(double cur_time, int k, priority_queue<JOB>& rq, double extra_time){
//     vector<JOB> cur_current_jobs(current_jobs.begin(), current_jobs.begin() + num_of_task);
//     vector<double> partitionDeadlines(num_of_task);
//     for(int i = 0; i < num_of_task; i++){
//         partitionDeadlines[i] = cur_current_jobs[i].deadline;
//     }
//     sort(cur_current_jobs.begin(), cur_current_jobs.end(), [](JOB a, JOB b) {
//         return a.deadline < b.deadline;
//     });
//     sort(partitionDeadlines.begin(), partitionDeadlines.end());

//     priority_queue<JOB> tmp;  // to store slackless-discarded jobs
//     double disc_used_slack = 0.0;
//     while(!discardedJobs.empty()){
//         JOB disc_job = discardedJobs.top();
//         discardedJobs.pop();
//         if(current_level > k){
//             disc_job.deadline = disc_job.original_deadline; // reseting deadline to original deadline
//         }
//         vector<vector<JOB>> jobPartitions(num_of_task);
//         for(auto job:jobs){
//             if(current_level > k){
//                 job.deadline = job.original_deadline; // resetting deadline to original deadline
//             }
//             if((job.deadline <= cur_time) or (job.criticality < current_level) or (job.arrival_time >= disc_job.deadline))
//                 continue;
//             int partition = findPartition(job, partitionDeadlines);
//             jobPartitions[partition].push_back(job);
//         }
//         double cur_slackSystem = INT_MAX;
//         for(int i = 0; i < num_of_task; i++){
//             int size = jobPartitions[i].size();
//             if(size == 0)
//                 continue;
//             // printf("Partition[%d] size = %d\n", i, size);
//             double cur_slack = slackTable[current_level-1][jobPartitions[i][0].deadline_index][jobPartitions[i][size-1].deadline_index] - total_idle_time - total_stolen_time - findSumExeCurrentJobs(i+1, cur_current_jobs);

//             // double cur_slack = findPreSlack(jobPartitions[i][0].deadline_index, jobPartitions[i][size-1].deadline_index, disc_job.deadline) - total_idle_time - total_stolen_time - findSumExeCurrentJobs(i+1, cur_current_jobs);

//             // printf("pre-slack([%d][%d]): %lf\n", jobPartitions[i][0].deadline_index, jobPartitions[i][size-1].deadline_index, slackTable[current_level-1][jobPartitions[i][0].deadline_index][jobPartitions[i][size-1].deadline_index]);
//             // printf("total_idle_time: %lf\n", total_idle_time);
//             // printf("total_stolen_time: %lf\n", total_stolen_time);
//             // printf("findSumExeCurrentJobs(i+1, cur_current_jobs): %lf\n\n", findSumExeCurrentJobs(i+1, cur_current_jobs));
//             cur_slackSystem = min(cur_slackSystem, cur_slack);  // ^ re-calculate execution times of each index in slackTable???
//         }
//         if(cur_slackSystem == INT_MAX){
//             tmp.push(disc_job);
//         }
//         else{
//             cur_slackSystem += extra_time;
//             cur_slackSystem = round(cur_slackSystem * 1e6) / 1e6;
//             printf("\ncur_slackSystem: %lf\n", cur_slackSystem);
//             if((cur_slackSystem - disc_used_slack) >= disc_job.remain_time){ 
//                 printf("\nT%d(J%d) with remain time:%lf undiscarded due to slack of %lf\n", disc_job.task_id + 1, disc_job.job_id + 1, disc_job.remain_time ,cur_slackSystem - disc_used_slack);
//                 disc_used_slack += cur_slackSystem;
//                 rq.push(disc_job);
//             }
//             else{
//                 tmp.push(disc_job);
//             }
//         }
//     }
//     discardedJobs = tmp;
// }

double dpm(double cur_time, priority_queue<JOB> rq, int k){
    double miniSlack = INT_MAX;
    priority_queue<JOB> rqCopy = rq;
    vector<JOB> rqJobs;
    while(!rqCopy.empty()){
        JOB job = rqCopy.top();
        rqCopy.pop();
        rqJobs.push_back(job);
    }
    for(auto job:jobs){
        if(job.arrival_time <= cur_time or job.criticality < current_level)
            continue;
        if(current_level > k){
            job.deadline = job.original_deadline; // resetting deadline to original deadline
        }
        rqJobs.push_back(job);
    }
    for(int i = 0; i < rqJobs.size(); i++){ // QOS for already undiscarded jobs but not for those which will be discarded later and might have been undiscarded due to slack
        JOB job = rqJobs[i];
        if(current_level > k){ // DEBUG
            job.deadline = job.original_deadline; // resetting deadline to original deadline
        }
        double slack = job.deadline - cur_time;
        for(int j = 0; j < rqJobs.size(); j++){
            JOB other_job = rqJobs[j]; // we need to consider execution time of current job also since we are finding slack of system due to this job and later checking if slack > 0
            if(other_job.arrival_time >= job.deadline)
                continue;
            if(current_level > k){ // DEBUG
                other_job.deadline = other_job.original_deadline; // resetting deadline to original deadline
            }
            double exeFrac = min(1.0, (job.deadline - cur_time)/(other_job.deadline - cur_time)); // there would be some undiscarded job also in readyQueue
            slack -= exeFrac * (execution_time[other_job.task_id][min(current_level, other_job.criticality)-1] - other_job.executed_time);
            // printf("T%d(J%d):%lf, ", other_job.task_id+1, other_job.job_id+1, slack);
        }
        // 10 - 6 - 4
        miniSlack = min(miniSlack, slack);
        miniSlack = round(miniSlack * 1e6) / 1e6;
        // printf("miniSlack[T%d(J%d)]: %lf\n", job.task_id+1, job.job_id+1, miniSlack);
    }
    return miniSlack;
}
void computeCurrentSlack2(double cur_time, int k, priority_queue<JOB>& rq, JOB cur_job, double extra_time){
    priority_queue<JOB> tmp;  // to store slackless-discarded jobs
    priority_queue<JOB> rqCopy = rq;
    double disc_used_slack = 0.0;
    while(!discardedJobs.empty()){
        JOB disc_job = discardedJobs.top();
        discardedJobs.pop();
        if(disc_job.deadline <= cur_time){
            dynamic_discarded++;
            continue;
        }
        if(current_level > k){
            disc_job.deadline = disc_job.original_deadline; // reseting deadline to original deadline
        }
        double slack = disc_job.deadline - cur_time;
        printf("\ndisc_job.deadline: %lf\n", disc_job.deadline);


        if(current_level > k){
            cur_job.deadline = cur_job.original_deadline; // resetting deadline to original deadline
        }
        // double exeFrac = min(1.0, (disc_job.deadline - cur_time)/(cur_job.deadline - cur_time));
        // slack -= exeFrac * (execution_time[cur_job.task_id][min(current_level, cur_job.criticality)-1] - cur_job.executed_time);
        // printf("%T%d(J%d):%lf, ", cur_job.task_id+1, cur_job.job_id+1, cur_job.deadline);

        while(!rqCopy.empty()){
            JOB job = rqCopy.top();
            rqCopy.pop();
            double exeFrac = min(1.0, (disc_job.deadline - cur_time)/(job.deadline - cur_time)); // there would be some undiscarded job also in readyQueue
            slack -= exeFrac * (execution_time[job.task_id][min(current_level, job.criticality)-1] - job.executed_time);
            // printf("^T%d(J%d):%lf, ", job.task_id+1, job.job_id+1, job.deadline);
        }
        for(auto job:jobs){
            if(job.arrival_time <= cur_time or job.arrival_time >= disc_job.deadline or job.criticality < current_level)
                continue;
            if(current_level > k){
                job.deadline = job.original_deadline; // resetting deadline to original deadline
            }
            double exeFrac = min(1.0, (disc_job.deadline - job.arrival_time)/(job.deadline - job.arrival_time));
            slack -= exeFrac * execution_time[job.task_id][current_level-1];
            // printf("T%d(J%d):%lf, ", job.task_id+1, job.job_id+1, job.deadline);

        }
        printf("\nslack: %lf\n", slack);
        if((slack - disc_used_slack) >= (execution_time[disc_job.task_id][disc_job.criticality-1] - disc_job.executed_time)){
            printf("\nT%d(J%d) with remain time:%lf undiscarded due to slack of %lf\n\n", disc_job.task_id + 1, disc_job.job_id + 1, execution_time[disc_job.task_id][disc_job.criticality-1] - disc_job.executed_time ,slack);
            total_stolen_time += slack - disc_used_slack;
            disc_used_slack += slack;
            rq.push(disc_job);
        }
        else{
            tmp.push(disc_job);
        }
    }
    discardedJobs = tmp;

}
void pushRqIntoDiscardedJobs(priority_queue<JOB>& pqLO){ 
    priority_queue<JOB> tmp;
    while(!pqLO.empty()){
        JOB cur_job = pqLO.top();
        pqLO.pop();
        if(cur_job.criticality < current_level){
            cur_job.discarded = true;
            discardedJobs.push(cur_job);
        }
        else{
            tmp.push(cur_job);
        }
    }
    pqLO = tmp;
}

void doSlackWorks(double& cur_time, int& k, priority_queue<JOB>& pqLO, JOB cur_job, double extra_time){
    if(discardedJobs.empty())
        return;
    // computeCurrentSlack(cur_time, k, pqLO, extra_time);
    computeCurrentSlack2(cur_time, k, pqLO, cur_job, extra_time);
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
            printf("NO SOLUTION EXISTS\n");
            k = max_criticality;   ///////// CHANGE
            exit(1);
        }
        double x = findX(k);
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

	JOB cur_job;
	double cur_time = 0.0, prev_time = 0.0, extra_time = 0.0;
	int cur_job_index = -1, prev_job_index = -1; // no job in CPU

	double decision_point = findMinTask();
    if(decision_point != cur_time){
        total_idle_time += decision_point - cur_time;
	    printf("%lf -> %lf => IDLE\n", cur_time, decision_point);
    }
	cur_time = decision_point;
	bool arrival = false;

    while(cur_time < time){
        bool idle = cur_job_index == -1;
        if(cur_job_index != -1){
            cur_job.remain_time -= cur_time - prev_time;
            cur_job.executed_time += cur_time - prev_time;

            cur_job.remain_time = round(cur_job.remain_time * 1e6) / 1e6;
            cur_job.executed_time = round(cur_job.executed_time * 1e6) / 1e6;
            current_jobs[cur_job.task_id] = cur_job;


            if(abs(cur_job.remain_time) <= 0){  // remove completed job from CPU
                // printf("Job %d of Task %d completed\n", cur_job.job_id + 1, cur_job.task_id + 1);                
                if(!cur_job.discarded){  
                    extra_time += execution_time[cur_job.task_id][current_level-1] - cur_job.executed_time;
                    if(cur_job.executed_time < COMPLETION_LB * execution_time[cur_job.task_id][current_level-1]){
                        printf("T%d(J%d) completed early\n", cur_job.task_id + 1, cur_job.job_id + 1);
                        doSlackWorks(cur_time, k, pqLO, cur_job, extra_time);
                    }
                }
                if(pqLO.empty()){
                    double procrast_slack = dpm(cur_time, pqLO, k);
                    if(procrast_slack > PROCRAST_THRESHOLD){
                        decision_point = min(cur_time + procrast_slack, findMinTask());
                        total_idle_time += decision_point - cur_time;
                        printf("%lf -> %lf => PROCRAST IDLE\n", cur_time, decision_point);
                        usleep(SLEEP_TIME);
                        prev_time = cur_time;
                        cur_time = decision_point;
                        arrival = false;
                        cur_job_index = -1; 
				        prev_job_index = -1;
                        continue;
                    }
                }
                cur_job_index = -1; 
				prev_job_index = -1;
			}

            if(!cur_job.discarded && (cur_job.executed_time >= execution_time[cur_job.task_id][current_level-1])){ // CHANGE >= to >
                printf("\nCriticality Level changed from %d to %d\n\n", current_level, current_level+1);
                current_level++;
                printf("Job %d of Task %d preempted(%lf)\n", cur_job.job_id + 1, cur_job.task_id + 1, cur_job.remain_time);
                cur_job.preempted_time = cur_time;  // current_level changing is also a preemption point right now
                pqLO.push(cur_job);
                cur_job_index = -1;
                pushRqIntoDiscardedJobs(pqLO);
                if(current_level > k){
                    transferRq(pqLO);
                }
                doSlackWorks(cur_time, k, pqLO, cur_job, extra_time);
            }
            
            if(cur_job_index != -1){
                // printf("Job %d of Task %d preempted(%lf)\n", cur_job.job_id + 1, cur_job.task_id + 1, cur_job.remain_time);
                cur_job.preempted_time = cur_time;  // current_level changing is also a preemption point right now
                pqLO.push(cur_job);
                cur_job_index = -1;
            }
        }
        if(!arrival){ 
			int discardedCount = handleArrivals(cur_time, k, pqLO);
            if(discardedCount > 0)
                doSlackWorks(cur_time, k, pqLO, cur_job, extra_time);
            if(idle){
                double procrast_slack = dpm(cur_time, pqLO, k);
                if(procrast_slack > PROCRAST_THRESHOLD){
                    decision_point = min(cur_time + procrast_slack, findMinTask());
                    total_idle_time += decision_point - cur_time;
                    printf("%lf -> %lf => PROCRAST IDLE\n", cur_time, decision_point);
                    usleep(SLEEP_TIME);
                    prev_time = cur_time;
                    cur_time = decision_point;
                    arrival = false;
                    continue;
                }

            }
			arrival = true;  
		}
        
		
		
		if(!pqLO.empty()){
            cur_job = pqLO.top();
            pqLO.pop();

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

            if(!cur_job.discarded)
			    decision_point = min({cur_time + cur_job.remain_time, cur_time + execution_time[cur_job.task_id][current_level-1] - cur_job.executed_time, findMinTask()});
            else    
                decision_point = min({cur_time + cur_job.remain_time, findMinTask()});

            
			if(decision_point == cur_time + cur_job.remain_time){
                if(!cur_job.discarded)
				    printf("%lf -> %lf => T%d(J%d)*\n", cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                else{
                    printf("%lf -> %lf => Discarded T%d(J%d)*\n", cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                    // printf("cur_job.remain_time: %lf\n", cur_job.remain_time);
                    // printf("execution_time[cur_job.task_id][current_level-1] = %lf\n", execution_time[cur_job.task_id][current_level-1]);
                }
			}
			else{
				if(!cur_job.discarded)
				    printf("%lf -> %lf => T%d(J%d)\n", cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                else{
                    printf("%lf -> %lf => Discarded T%d(J%d)\n", cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                    // printf("cur_job.remain_time: %lf\n", cur_job.remain_time);
                    // printf("execution_time[cur_job.task_id][current_level-1] = %lf\n", execution_time[cur_job.task_id][current_level-1]);
                }
			}
		}
		else{
			decision_point = findMinTask();
            total_idle_time += decision_point - cur_time;
			printf("%lf -> %lf => IDLE\n", cur_time, decision_point);
		}
        usleep(SLEEP_TIME);
		prev_time = cur_time;
		cur_time = decision_point;
		arrival = false;

    }
}

void printSlackTable() {
    for (int k = 0; k < max_criticality; k++) {
        for (int i = 0; i < num_of_jobs; i++) {
            for (int j = 0; j < num_of_jobs; j++) {
                printf("%.1lf ", slackTable[k][i][j]);
            }
            printf("\n");
        }
        printf("-----------------------------------\n");
    }
}

void printJobs() {
    for (auto job:jobs) {
        printf("T%d(J%d) : Deadline = %lf\n", job.task_id+1, job.job_id+1, job.deadline);
    }
}
  
void edfvd(int hyperperiod) {
    int k = offlinePreprocessing();
    printf("k: %d\n", k);
    loadJobsInHyperperiod(hyperperiod);
    // printJobs();
    num_of_jobs = jobs.size();
    slackTablePreCompute();
    // printSlackTable();
    runtimeScheduling(k, hyperperiod);
}

int main(int argc, char* argv[]){
    get_task_info();		//collecting tasks detail
	printf("Hyperperiod: %d\n", hyperperiod);
    edfvd(hyperperiod);
	printf("Preemptions: %d\n", preemptions);
    printf("Discarded Jobs: %d\n", dynamic_discarded);
    printf("Current Discarded Queue Size: %d\n", discardedJobs.size());
	return 0;
}





/*
2
0 4 4 1 2
0 6 6 2 1 5
*/

/* 
3
0 10 10 1 4
0 12 12 2 2 4
0 15 15 3 1 3 4
*/

/*
BEST
3
0 12 12 1 4
0 20 20 2 2 4
0 10 10 3 3 5 6
*/

// What if a discarded job has slack of 5 at time t = 10 and current_critical_level, executes the remaining execution time and exits but then the HI criticality job changes criticality. Won't the previously calculated slack for discarded job based on lower critical level wcet of HI criticality job be wrong and create negative slack now? Or is edf-vd taking care of it?