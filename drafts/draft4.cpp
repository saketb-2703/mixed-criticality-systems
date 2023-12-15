#include <bits/stdc++.h>
#include<unistd.h>  
using namespace std;


/*
The logic is implemented using 3 types of threads depecting the 3 types of cores:
1. Non-shutdownable cores
2. Shutdownable cores
3. Exceptional core

The Non-shutdownable cores use only DVFS to conserve energy. In the rare case the utilisation happens to go beyond a certain point, it will try to
migrate jobs to other non-shutdownable cores.

The Shutdownable used DPM and DVFS to conserve energy. 
*/
#define MAX_TASK 10
#define MAX_JOBS 100
#define NO_CORES 5
#define MAX_CRITICALITY 10
#define SLEEP_TIME 10000
#define COMPLETION_LB 0.8
#define PROCRAST_THRESHOLD 0.5
struct job{
	double arrival_time;      // clock time at which this job was admitted in the system
	int period;              // Period instantiated from the task set
	double executed_time;   // So far executed time
    double executed_time_hfq; // So far executed time if job had been executed at highest frequency throughout
	double deadline = -1;     // Absolute deadline of the job(will be virtual absolute deadline if current criticality <= k)
    double original_deadline; // Absolute deadline of the job(ignoring EDF-VD virtual deadline)
	int task_id = -1;         // task_id of job from task set
	int job_id = -1;          // job_index of job considering this task only
	int job_index;            // job_index of job considering all jobs arrived so far of all tasks
    int deadline_index;       // job_index of job when all the jobs are sorted by their absolute deadlines
	double remain_time;       // remaining actual execution time of this job; instantiated at runtime when the job is admitted in the system
	double relative_deadline; // relative deadline of the job as per task set
	double cur_time;          // clock time at which this job was admitted in the system   // same as arrival time..so ignore
	double preempted_time = -1; // clock time at which this job was preempted from the core
    int criticality;           // criticality level of this job(i.e. task)
    bool discarded = false;    // indicated whether this job was part of discarded queue and hence later admitted into the ready queue 

	bool operator<(const job& other) const {
		if(deadline == other.deadline) 
			return preempted_time < other.preempted_time; // The latest preempted/cpu-executed job should be selected among equal priority jobs
        return deadline > other.deadline; 
    }
};

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
int global_num_of_task;
typedef struct job JOB;
vector<custom_priority_queue<JOB>> rq(NO_CORES); //vector of ready queues in which i-th element is the ready queue for i-th core
double global_arrival_time[MAX_TASK], global_execution_time[MAX_TASK][MAX_CRITICALITY], global_deadline[MAX_TASK], global_virtualDeadlines[MAX_TASK], global_slackTable[MAX_CRITICALITY][MAX_JOBS][MAX_JOBS];
int global_criticality[MAX_TASK], global_period[MAX_TASK];
//data structures assosciated with each core for each criticality level and global data structures for all cores
vector<bool> shutdownable(NO_CORES);  //boolean array initialised in the beginning to indicate whether a core is shutdownable or not
vector<vector<JOB> > migratedJobs(NO_CORES); /*global data structure which indicates the jobs migrated from a core to another core
the i-th element of this vector is a vector of jobs migrated from other cores to the t-th core. The i-th core has to accomodate all these jobs into its ready queue*/    
vector<JOB> cur_jobs(NO_CORES); // stores current executing job of each core
vector<int> current_level(NO_CORES, 1); // stores current criticality level of each core

class NonShutDownableCore
{
    public:
        int coreNumber; //coreNumber of current core
        int num_of_task;
        int num_of_jobs;
        int hyperperiod;
        int job_index = 0;
        int preemptions = 0;
        int max_criticality = 0;
        int dynamic_discarded = 0;
        vector<double> utilTask{vector<double>(MAX_TASK,0.0)};  //todo
        vector<double> availableFrequency = {0.5, 0.75, 1.0};  //hardcoded array with the available frequencies of the core
        double Foptimal_min = 0.55, Foptimal_max = 0.85; //hardcoded optimal frequency range for the core
        vector<JOB> jobs;
        double total_idle_time = 0.0, total_stolen_time = 0.0, slackSystem = 0.0; //todo
        int cur_job_index = -1, prev_job_index = -1; // no job in CPU
        vector<JOB> current_jobs{vector<JOB>(MAX_TASK)};
        custom_priority_queue<JOB> discardedJobs, pendingJobs;
        double arrival_time[MAX_TASK], execution_time[MAX_TASK][MAX_CRITICALITY], deadline[MAX_TASK], virtualDeadlines[MAX_TASK], slackTable[MAX_CRITICALITY][MAX_JOBS][MAX_JOBS];
        int criticality[MAX_TASK], period[MAX_TASK], global_id[MAX_TASK]; //data structures to store task set information of this core
        vector<int> cur_task_info{vector<int>(MAX_TASK, 0)}; 

        void mainEq(vector<int> myTasks, int coreNo) //entry point in thread
        {
            cout << "[Core " << coreNo << "]" << endl;
            coreNumber = coreNo;
            num_of_task = myTasks.size();
            init(myTasks);
            edfvd(hyperperiod); //actual implementation starts
            printf("Preemptions: %d\n", preemptions);
            printf("Discarded Jobs: %d\n", dynamic_discarded);
            printf("Current Discarded Queue Size: %ld\n", discardedJobs.size());


        }
        JOB loadJob(int task_id, double cur_time){  //create a job object for the given task, initialise the required data and return the object
            JOB t;
            t.arrival_time = cur_time;
            t.period = period[task_id];
            t.executed_time = 0;
            t.executed_time_hfq = 0;
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
        double findMinTask(){  //find the next job arrival
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
        bool checkKCondition(int k){  //todo
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
        int findK(){ //todo
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
        double findX(int k){ //todo
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

        void transferRq(custom_priority_queue<JOB>& pqLO){ 
            custom_priority_queue<JOB> tmp;
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


        // void slackTablePreCompute(){
        //     vector<vector<double>> initialSlack(max_criticality, vector<double>(num_of_jobs));
        //     for(int p = 1; p <= max_criticality; p++){
        //         for (int i = 0; i < num_of_jobs; i++) {
        //             initialSlack[p-1][i] = jobs[i].deadline - 0; // Assuming time origin is 0
        //             for (int j = 0; j < num_of_jobs; j++) {
        //                 if((jobs[j].criticality >= p) && (jobs[j].deadline <= jobs[i].deadline))  // need to consider execution times of jobs with criticality >= p[because other jobs are discarded] and deadline less than or equal to deadline of job i
        //                     initialSlack[p-1][i] -= execution_time[jobs[j].task_id][p-1];
        //             }
        //             // cout << initialSlack[p-1][i] << ",";
        //         }
        //         // cout << endl;
        //     }
        //     for(int p = 1; p <= max_criticality; p++){
        //         for(int i = 0; i < num_of_jobs; i++){
        //             double miniSlack = INT_MAX;
        //             for(int j = i; j < num_of_jobs; j++){
        //                 miniSlack = min(miniSlack, initialSlack[p-1][j]);
        //                 slackTable[p-1][i][j] = miniSlack;
        //             }
        //         }
        //     }
        // }

        int handleArrivals(double cur_time, int k, custom_priority_queue<JOB>& pqLO){
            int HIcount = 0, LOcount = 0;
            for (int j = 0; j < num_of_task; j++) {
                if (fmod((cur_time - arrival_time[j]), period[j]) == 0) {
                    JOB t = loadJob(j, cur_time);
                    if(t.criticality >= current_level[coreNumber]){
                        pqLO.push(t);
                        HIcount++;
                    }
                    else{
                        t.discarded = true;
                        pendingJobs.push(t);
                        LOcount++;
                    }
                    // printf("Job %d of Task %d arrived with criticality %d\n", t.job_id + 1, t.task_id + 1, t.criticality);
                    cur_task_info[j]++;
                }
            }
            return HIcount;
            
        }

        double dpm(double cur_time, custom_priority_queue<JOB> rq, int k){
            double miniSlack = INT_MAX;
            custom_priority_queue<JOB> rqCopy = rq;
            vector<JOB> rqJobs;
            while(!rqCopy.empty()){
                JOB job = rqCopy.top();
                rqCopy.pop();
                rqJobs.push_back(job);
            }
            for(auto job:jobs){
                if(job.arrival_time <= cur_time or job.criticality < current_level[coreNumber])
                    continue;
                if(current_level[coreNumber] > k){
                    job.deadline = job.original_deadline; // resetting deadline to original deadline
                }
                rqJobs.push_back(job);
            }
            for(int i = 0; i < rqJobs.size(); i++){ // QOS for already undiscarded jobs but not for those which will be discarded later and might have been undiscarded due to slack
                JOB job = rqJobs[i];
                if(current_level[coreNumber] > k){ // DEBUG
                    job.deadline = job.original_deadline; // resetting deadline to original deadline
                }
                double slack = job.deadline - cur_time;
                for(int j = 0; j < rqJobs.size(); j++){
                    JOB other_job = rqJobs[j]; // we need to consider execution time of current job also since we are finding slack of system due to this job and later checking if slack > 0
                    if(other_job.arrival_time >= job.deadline)
                        continue;
                    if(current_level[coreNumber] > k){ // DEBUG
                        other_job.deadline = other_job.original_deadline; // resetting deadline to original deadline
                    }
                    double exeFrac = min(1.0, (job.deadline - cur_time)/(other_job.deadline - cur_time)); // there would be some undiscarded job also in readyQueue
                    slack -= exeFrac * (execution_time[other_job.task_id][min(current_level[coreNumber], other_job.criticality)-1] - other_job.executed_time_hfq);
                    // printf("T%d(J%d):%lf, ", other_job.task_id+1, other_job.job_id+1, slack);
                }
                // 10 - 6 - 4
                miniSlack = min(miniSlack, slack);
                miniSlack = round(miniSlack * 1e6) / 1e6;
                // printf("miniSlack[T%d(J%d)]: %lf\n", job.task_id+1, job.job_id+1, miniSlack);
            }
            return miniSlack;
        }
        void computeCurrentSlack2(double cur_time, int k, custom_priority_queue<JOB>& rq, JOB cur_job, double extra_time){ //func to calculate slack
            custom_priority_queue<JOB> tmp;  // to store slackless-discarded jobs
            custom_priority_queue<JOB> rqCopy = rq;
            double disc_used_slack = 0.0;
            while(!discardedJobs.empty()){
                JOB disc_job = discardedJobs.top();
                discardedJobs.pop();
                if(disc_job.deadline <= cur_time){
                    dynamic_discarded++;
                    continue;
                }
                if(current_level[coreNumber] > k){
                    disc_job.deadline = disc_job.original_deadline; // reseting deadline to original deadline
                }
                double slack = disc_job.deadline - cur_time;
                // printf("\ndisc_job.deadline: %lf\n", disc_job.deadline);


                if(current_level[coreNumber] > k){
                    cur_job.deadline = cur_job.original_deadline; // resetting deadline to original deadline
                }
                // double exeFrac = min(1.0, (disc_job.deadline - cur_time)/(cur_job.deadline - cur_time));
                // slack -= exeFrac * (execution_time[cur_job.task_id][min(current_level[coreNumber], cur_job.criticality)-1] - cur_job.executed_time);
                // printf("%T%d(J%d):%lf, ", cur_job.task_id+1, cur_job.job_id+1, cur_job.deadline);

                while(!rqCopy.empty()){
                    JOB job = rqCopy.top();
                    rqCopy.pop();
                    double exeFrac = min(1.0, (disc_job.deadline - cur_time)/(job.deadline - cur_time)); // there would be some undiscarded job also in readyQueue
                    slack -= exeFrac * (execution_time[job.task_id][min(current_level[coreNumber], job.criticality)-1] - job.executed_time_hfq);
                    // printf("^T%d(J%d):%lf, ", job.task_id+1, job.job_id+1, job.deadline);
                }
                for(auto job:jobs){
                    if(job.arrival_time <= cur_time or job.arrival_time >= disc_job.deadline or job.criticality < current_level[coreNumber])
                        continue;
                    if(current_level[coreNumber] > k){
                        job.deadline = job.original_deadline; // resetting deadline to original deadline
                    }
                    double exeFrac = min(1.0, (disc_job.deadline - job.arrival_time)/(job.deadline - job.arrival_time));
                    slack -= exeFrac * execution_time[job.task_id][current_level[coreNumber]-1];
                    // printf("T%d(J%d):%lf, ", job.task_id+1, job.job_id+1, job.deadline);

                }
                printf("\n[Core %d] slack: %lf\n", coreNumber, slack);
                if((slack - disc_used_slack) >= (execution_time[disc_job.task_id][disc_job.criticality-1] - disc_job.executed_time_hfq)){
                    printf("\n[Core %d] T%d(J%d) with remain time:%lf undiscarded due to slack of %lf\n\n", coreNumber, disc_job.task_id + 1, disc_job.job_id + 1, execution_time[disc_job.task_id][disc_job.criticality-1] - disc_job.executed_time_hfq ,slack);
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
        void pushRqIntoDiscardedJobs(custom_priority_queue<JOB>& pqLO){ //
            custom_priority_queue<JOB> tmp;
            while(!pqLO.empty()){
                JOB cur_job = pqLO.top();
                pqLO.pop();
                if(cur_job.criticality < current_level[coreNumber]){
                    cur_job.discarded = true;
                    discardedJobs.push(cur_job);
                }
                else{
                    tmp.push(cur_job);
                }
            }
            pqLO = tmp;
        }

        void doSlackWorks(double& cur_time, int& k, custom_priority_queue<JOB>& pqLO, JOB cur_job, double extra_time){
            if(pendingJobs.empty())
                return;
            while(!pendingJobs.empty()){
                JOB job = pendingJobs.top();
                pendingJobs.pop();
                discardedJobs.push(job);
            }
            // computeCurrentSlack(cur_time, k, pqLO, extra_time);
            computeCurrentSlack2(cur_time, k, pqLO, cur_job, extra_time);
        }

        double calculateFrequency(){  //DVFS algo for frequency selection
            double sum = 0.0;
            for(int i = 0; i < num_of_task; i++){
                sum += utilTask[i];
            }
            // printf("sum = %lf\n", sum);
            int index = upper_bound(availableFrequency.begin(), availableFrequency.end(), sum) - availableFrequency.begin();
            return availableFrequency[min(index, (int)availableFrequency.size()-1)];
        }

        double computeCoreSlack(double cur_time, int k, int core, custom_priority_queue<JOB> rq, JOB cur_job, JOB avoid_job){
            double miniSlack = INT_MAX;
            custom_priority_queue<JOB> rqCopy = rq;
            if(cur_job.job_index != avoid_job.job_index)
                rqCopy.push(cur_job);
            vector<JOB> rqJobs;
            while(!rqCopy.empty()){
                JOB job = rqCopy.top();
                rqCopy.pop();
                rqJobs.push_back(job);
            }
            for(auto job:jobs){
                if(job.arrival_time <= cur_time or job.criticality < current_level[core] or job.job_index == avoid_job.job_index)
                    continue;
                if(current_level[core] > k){
                    job.deadline = job.original_deadline; // resetting deadline to original deadline
                }
                rqJobs.push_back(job);
            }
            for(int i = 0; i < rqJobs.size(); i++){ // QOS for already undiscarded jobs but not for those which will be discarded later and might have been undiscarded due to slack
                JOB job = rqJobs[i];
                if(current_level[core] > k){ // DEBUG
                    job.deadline = job.original_deadline; // resetting deadline to original deadline
                }
                double slack = job.deadline - cur_time;
                for(int j = 0; j < rqJobs.size(); j++){
                    JOB other_job = rqJobs[j]; // we need to consider execution time of current job also since we are finding slack of system due to this job and later checking if slack > 0
                    if(other_job.arrival_time >= job.deadline)
                        continue;
                    if(current_level[core] > k){ // DEBUG
                        other_job.deadline = other_job.original_deadline; // resetting deadline to original deadline
                    }
                    double exeFrac = min(1.0, (job.deadline - cur_time)/(other_job.deadline - cur_time)); // there would be some undiscarded job also in readyQueue
                    slack -= exeFrac * (execution_time[other_job.task_id][min(current_level[core], other_job.criticality)-1] - other_job.executed_time_hfq);
                    // printf("T%d(J%d):%lf, ", other_job.task_id+1, other_job.job_id+1, slack);
                }
                // 10 - 6 - 4
                miniSlack = min(miniSlack, slack);
                miniSlack = round(miniSlack * 1e6) / 1e6;
                // printf("miniSlack[T%d(J%d)]: %lf\n", job.task_id+1, job.job_id+1, miniSlack);
            }
            return miniSlack;
        }
        bool migratable(double cur_time, int migratedCore, JOB mig_job, int k){
            custom_priority_queue<JOB> coreRqCopy = rq[migratedCore];
            coreRqCopy.push(mig_job); // assuming migrating job inside ready queue of migrated core(eventhough it didn't arrive yet in real time)
            JOB dummy_job;
            dummy_job.job_index = -1;
            double slack = computeCoreSlack(cur_time, k, migratedCore, coreRqCopy, cur_jobs[migratedCore], dummy_job); // if u call doSlackWorks, then pendingJobsQueue empty is also considered
            if(slack > 0)
                return true;
            return false;
        }

        void migrateJobs(double cur_time, int k){ //migrate low criticality jobs from this core to available core
            ////////////////////////////////////////////// migrate either upcoming job of every task or jobs in ready queue
            double maxSlack = INT_MIN;
            JOB jobToRemove;
            int typeJobToRemove = -1; // 0 - RQ job, 1 - upcoming job, 2 - current executing job
            custom_priority_queue<JOB> rqCopy = rq[coreNumber];
            vector<JOB> rqJobs;
            while(!rqCopy.empty()){
                JOB job = rqCopy.top();
                rqCopy.pop();
                rqJobs.push_back(job);
            }
            rqCopy = rq[coreNumber];
            for(auto job:rqJobs){
                rqCopy.remove(job);
                JOB dummy_job;
                dummy_job.job_index = -1;
                double slack = computeCoreSlack(cur_time, k, coreNumber, rqCopy, cur_jobs[coreNumber], dummy_job); // slack obtained by removing RQ job
                if (slack > maxSlack) {
                    maxSlack = slack;
                    jobToRemove = job;
                    typeJobToRemove = 0;
                }
                rqCopy.push(job); // Put back the previously removed job
            }

            int num_of_upcoming_checks = num_of_task;
            for(auto job:jobs){
                if(job.arrival_time <= cur_time)
                    continue;
                double slack = computeCoreSlack(cur_time, k, coreNumber, rqCopy, cur_jobs[coreNumber], job); // slack obtained by removing upcoming job
                if (slack > maxSlack) {
                    maxSlack = slack;
                    jobToRemove = job;
                    typeJobToRemove = 1;
                }
                num_of_upcoming_checks--;
                if(num_of_upcoming_checks == 0)
                    break;
            }
            double slack = computeCoreSlack(cur_time, k, coreNumber, rqCopy, cur_jobs[coreNumber], cur_jobs[coreNumber]); // slack obtained by removing current executing job in core
            if (slack > maxSlack) {
                maxSlack = slack;
                jobToRemove = cur_jobs[coreNumber];
                typeJobToRemove = 2;
            }

            ////////////////////////////////////////////////
            JOB mig_job;
            int task_id_jobToRemove = jobToRemove.task_id;
            if(typeJobToRemove == 0){
                mig_job = jobToRemove;
                rq[coreNumber].remove(jobToRemove);
            }
            else if(typeJobToRemove == 1){
                mig_job = loadJob(task_id_jobToRemove, (cur_task_info[task_id_jobToRemove]) * period[task_id_jobToRemove] + arrival_time[task_id_jobToRemove]); // cur_time = (cur_task_info[i]) * period[i] + arrival_time[i]
            }
            else{
                mig_job = jobToRemove;
                cur_job_index = -1;
            }
            
            bool flag = false;
            for(int i = 0; i < NO_CORES; i++){
                if(i == coreNumber or shutdownable[i] or !migratable(cur_time, i, mig_job, k)) // check for LO criticality also
                    continue;
                migratedJobs[i].push_back(mig_job);
                flag = true;
                break;
            }
            if(flag){
                cur_task_info[task_id_jobToRemove]++; // this job has arrived in migrated core at a later stage, hence update as arrived already in current core.
                migrateJobs(cur_time , k); //recursive call to migrate other jobs if we were able to migrate this job
            }
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

        void runtimeScheduling(int k, int time){ //EDF-VD algo
            //custom_priority_queue<JOB> pqLO; // Min heap of pairs, where first is deadline and second is task id

            JOB cur_job;
            double cur_time = 0.0, prev_time = 0.0, extra_time = 0.0, freq = 1.0;
            

            double decision_point = findMinTask(); // find next job arrival
            if(decision_point != cur_time){
                total_idle_time += decision_point - cur_time;
                printf("[Core %d] %lf -> %lf => IDLE\n", coreNumber,cur_time, decision_point);
            }
            cur_time = decision_point;
            bool arrival = false;

            while(cur_time < time){
                bool idle = cur_job_index == -1;
                if(cur_job_index != -1){
                    cur_job.remain_time -= (cur_time - prev_time) * freq;
                    cur_job.executed_time_hfq += (cur_time - prev_time) * freq;
                    cur_job.executed_time += cur_time - prev_time;

                    cur_job.remain_time = round(cur_job.remain_time * 1e6) / 1e6;
                    cur_job.executed_time_hfq = round(cur_job.executed_time_hfq * 1e6) / 1e6;
                    cur_job.executed_time = round(cur_job.executed_time * 1e6) / 1e6;
                    current_jobs[cur_job.task_id] = cur_job;


                    if(abs(cur_job.remain_time) <= 0){  // remove completed job from CPU
                        // printf("Job %d of Task %d completed\n", cur_job.job_id + 1, cur_job.task_id + 1);                
                        utilTask[cur_job.task_id] = cur_job.executed_time_hfq / cur_job.deadline;
                        if(!cur_job.discarded){  


                            extra_time += execution_time[cur_job.task_id][current_level[coreNumber]-1] - cur_job.executed_time_hfq;
                            if(cur_job.executed_time < COMPLETION_LB * execution_time[cur_job.task_id][current_level[coreNumber]-1]){   // DOUBT
                                printf("[Core %d] T%d(J%d) completed early\n", coreNumber, cur_job.task_id + 1, cur_job.job_id + 1);
                                doSlackWorks(cur_time, k, rq[coreNumber], cur_job, extra_time);
                            }
                        }
                        if(rq[coreNumber].empty()){
                            migrateJobs(cur_time, k);
                            double procrast_slack = dpm(cur_time, rq[coreNumber], k);
                            if(procrast_slack > PROCRAST_THRESHOLD){
                                decision_point = min(cur_time + procrast_slack, findMinTask());
                                total_idle_time += decision_point - cur_time;
                                printf("[Core %d] %lf -> %lf => PROCRAST IDLE\n", coreNumber, cur_time, decision_point);
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

                    if(!cur_job.discarded && (cur_job.executed_time_hfq >= execution_time[cur_job.task_id][current_level[coreNumber]-1])){ // CHANGE >= to >
                        printf("\n[Core %d] Criticality Level changed from %d to %d\n\n", coreNumber, current_level[coreNumber], current_level[coreNumber]+1);
                        current_level[coreNumber]++;
                        printf("[Core %d]Job %d of Task %d preempted(%lf)\n", coreNumber, cur_job.job_id + 1, cur_job.task_id + 1, cur_job.remain_time);
                        cur_job.preempted_time = cur_time;  // current_level[coreNumber] changing is also a preemption point right now
                        rq[coreNumber].push(cur_job);
                        cur_job_index = -1;
                        pushRqIntoDiscardedJobs(rq[coreNumber]);
                        if(current_level[coreNumber] > k){
                            transferRq(rq[coreNumber]);
                        }
                        doSlackWorks(cur_time, k, rq[coreNumber], cur_job, extra_time);
                    }
                    
                    if(cur_job_index != -1){
                        // printf("Job %d of Task %d preempted(%lf)\n", cur_job.job_id + 1, cur_job.task_id + 1, cur_job.remain_time);
                        cur_job.preempted_time = cur_time;  // current_level[coreNumber] changing is also a preemption point right now
                        rq[coreNumber].push(cur_job);
                        cur_job_index = -1;
                    }
                }
                if(!arrival){ 
                    int HIcount = handleArrivals(cur_time, k, rq[coreNumber]);
                    if(HIcount > 0)
                        doSlackWorks(cur_time, k, rq[coreNumber], cur_job, extra_time);
                    if(idle){
                        migrateJobs(cur_time, k);  // what if a new job arrives when core is shutdown; it could be migrated to other core
                        double procrast_slack = dpm(cur_time, rq[coreNumber], k);
                        if(procrast_slack > PROCRAST_THRESHOLD){
                            decision_point = min(cur_time + procrast_slack, findMinTask());
                            total_idle_time += decision_point - cur_time;
                            printf("[Core %d] %lf -> %lf => PROCRAST IDLE\n", coreNumber, cur_time, decision_point);
                            usleep(SLEEP_TIME);
                            prev_time = cur_time;
                            cur_time = decision_point;
                            arrival = false;
                            continue;
                        }

                    }
                    arrival = true;  
                }
                
                
                
                if(!rq[coreNumber].empty()){
                    cur_job = rq[coreNumber].top();
                    rq[coreNumber].pop();

                    cur_job_index = cur_job.job_index;

                    if(cur_job.original_deadline <= cur_time){
                        printf("[Core %d] Job %d of Task %d missed deadline(%lf)\n", coreNumber, cur_job.job_id + 1, cur_job.task_id + 1, cur_job.deadline);
                        cur_job_index = -1;
                        continue;
                    }
                    if(cur_job.original_deadline - cur_time < cur_job.remain_time){
                        printf("[Core %d] Job %d of Task %d will miss deadline\n", coreNumber, cur_job.job_id + 1, cur_job.task_id + 1);
                        cur_job_index = -1;
                        continue;
                    }

                    if(prev_job_index != -1 && cur_job_index != prev_job_index){
                        preemptions++;
                        printf("[Core %d] Preemption at t = %lf\n", coreNumber, cur_time);
                    }
                    prev_job_index = cur_job_index;





                    utilTask[cur_job.task_id] = (1.0*execution_time[cur_job.task_id][criticality[cur_job.task_id]-1])/cur_job.deadline;
                    freq = calculateFrequency();
                    printf("\n[Core %d] Frequency changed to %lf\n", coreNumber, freq);
                    printf("[Core %d] Old remain time: %lf, ", coreNumber, cur_job.remain_time);
                    // cur_job.remain_time /= freq;
                    printf("[Core %d] New remain time: %lf\n", coreNumber, cur_job.remain_time / freq);

                    if(!cur_job.discarded)
                        decision_point = min({cur_time + cur_job.remain_time / freq, cur_time + execution_time[cur_job.task_id][current_level[coreNumber]-1] - cur_job.executed_time_hfq, findMinTask()});
                    else    
                        decision_point = min({cur_time + cur_job.remain_time / freq, findMinTask()});

                    
                    if(decision_point == cur_time + cur_job.remain_time){
                        if(!cur_job.discarded)
                            printf("[Core %d] %lf -> %lf => T%d(J%d)*\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                        else{
                            printf("[Core %d] %lf -> %lf => Discarded T%d(J%d)*\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                            // printf("cur_job.remain_time: %lf\n", cur_job.remain_time);
                            // printf("execution_time[cur_job.task_id][current_level[coreNumber]-1] = %lf\n", execution_time[cur_job.task_id][current_level[coreNumber]-1]);
                        }
                    }
                    else{
                        if(!cur_job.discarded)
                            printf("[Core %d] %lf -> %lf => T%d(J%d)\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                        else{
                            printf("[Core %d] %lf -> %lf => Discarded T%d(J%d)\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                            // printf("cur_job.remain_time: %lf\n", cur_job.remain_time);
                            // printf("execution_time[cur_job.task_id][current_level[coreNumber]-1] = %lf\n", execution_time[cur_job.task_id][current_level[coreNumber]-1]);
                        }
                    }
                }
                else{
                    decision_point = findMinTask();
                    total_idle_time += decision_point - cur_time;
                    printf("[Core %d] %lf -> %lf => IDLE\n", coreNumber, cur_time, decision_point);
                }
                usleep(SLEEP_TIME);
                prev_time = cur_time;
                cur_time = decision_point;
                arrival = false;

            }
        }

        void initUtilTask(){ //initilialise the utiltask structure
            for(int i = 0; i < num_of_task; i++){
                utilTask[i] = (1.0*execution_time[i][criticality[i]-1])/deadline[i];
                printf("[Core %d] utilTask[%d] = %lf\n", coreNumber, i, utilTask[i]);
            }
        }
        void printSlackTable() {  //print slack table
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
                printf("[Core %d] T%d(J%d) : Deadline = %lf\n", coreNumber, job.task_id+1, job.job_id+1, job.deadline);
            }
        }
        void edfvd(int hyperperiod) {
            int k = offlinePreprocessing(); //todo
            printf("k: %d\n", k);
            loadJobsInHyperperiod(hyperperiod); //todo
            // printJobs();
            num_of_jobs = jobs.size();
            // slackTablePreCompute();
            // printSlackTable();
            initUtilTask(); //initialise utilcal array which is used to calculate operating frequency during DVFS
            runtimeScheduling(k, hyperperiod); //start the execution of the tasks
        }
        void init(vector<int> myTasks)  //update the task information in the core from global information
        {
            for(int i = 0; i < num_of_task; i++)
            {
                global_id[i] = myTasks[i];
                arrival_time[i] = global_arrival_time[myTasks[i]];
                deadline[i] = global_deadline[myTasks[i]];
                criticality[i] = global_criticality[myTasks[i]];
                period[i] = global_period[myTasks[i]];
                virtualDeadlines[i] = global_virtualDeadlines[myTasks[i]];
                max_criticality = max(max_criticality, criticality[i]);
                for(int j = 0; j < criticality[i]; j++)
                {
                    execution_time[i][j] = global_execution_time[myTasks[i]][j];
                }
            }
            hyperperiod = findLCM();
        }
};
 
 // complete this thread.
class ShutDownableCore
{
    public:
        int coreNumber;
        int num_of_task;
        int num_of_jobs;
        int hyperperiod;
        int job_index = 0;
        int preemptions = 0;
        int max_criticality = 0;
        int dynamic_discarded = 0;
        vector<double> utilTask{vector<double>(MAX_TASK,0.0)};
        vector<double> availableFrequency = {0.5, 0.75, 1.0};
        double Foptimal_min = 0.55, Foptimal_max = 0.85;
        vector<JOB> jobs;
        double total_idle_time = 0.0, total_stolen_time = 0.0, slackSystem = 0.0;
        int cur_job_index = -1, prev_job_index = -1; // no job in CPU
        vector<JOB> current_jobs{vector<JOB>(MAX_TASK)};
        custom_priority_queue<JOB> discardedJobs, pendingJobs;
        double arrival_time[MAX_TASK], execution_time[MAX_TASK][MAX_CRITICALITY], deadline[MAX_TASK], virtualDeadlines[MAX_TASK], slackTable[MAX_CRITICALITY][MAX_JOBS][MAX_JOBS];
        int criticality[MAX_TASK], period[MAX_TASK], global_id[MAX_TASK];
        vector<int> cur_task_info{vector<int>(MAX_TASK, 0)};

        void mainEq(vector<int> myTasks, int coreNo)
        {
            cout << "[Core " << coreNo << "]" << endl;
            coreNumber = coreNo;
            num_of_task = myTasks.size();
            init(myTasks);
            edfvd(hyperperiod);
            printf("Preemptions: %d\n", preemptions);
            printf("Discarded Jobs: %d\n", dynamic_discarded);
            printf("Current Discarded Queue Size: %ld\n", discardedJobs.size());


        }
        JOB loadJob(int task_id, double cur_time){
            JOB t;
            t.arrival_time = cur_time;
            t.period = period[task_id];
            t.executed_time = 0;
            t.executed_time_hfq = 0;
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

        void transferRq(custom_priority_queue<JOB>& pqLO){
            custom_priority_queue<JOB> tmp;
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


        // void slackTablePreCompute(){
        //     vector<vector<double>> initialSlack(max_criticality, vector<double>(num_of_jobs));
        //     for(int p = 1; p <= max_criticality; p++){
        //         for (int i = 0; i < num_of_jobs; i++) {
        //             initialSlack[p-1][i] = jobs[i].deadline - 0; // Assuming time origin is 0
        //             for (int j = 0; j < num_of_jobs; j++) {
        //                 if((jobs[j].criticality >= p) && (jobs[j].deadline <= jobs[i].deadline))  // need to consider execution times of jobs with criticality >= p[because other jobs are discarded] and deadline less than or equal to deadline of job i
        //                     initialSlack[p-1][i] -= execution_time[jobs[j].task_id][p-1];
        //             }
        //             // cout << initialSlack[p-1][i] << ",";
        //         }
        //         // cout << endl;
        //     }
        //     for(int p = 1; p <= max_criticality; p++){
        //         for(int i = 0; i < num_of_jobs; i++){
        //             double miniSlack = INT_MAX;
        //             for(int j = i; j < num_of_jobs; j++){
        //                 miniSlack = min(miniSlack, initialSlack[p-1][j]);
        //                 slackTable[p-1][i][j] = miniSlack;
        //             }
        //         }
        //     }
        // }

        int handleArrivals(double cur_time, int k, custom_priority_queue<JOB>& pqLO){
            int HIcount = 0, LOcount = 0;
            for (int j = 0; j < num_of_task; j++) {
                if (fmod((cur_time - arrival_time[j]), period[j]) == 0) {
                    JOB t = loadJob(j, cur_time);
                    if(t.criticality >= current_level[coreNumber]){
                        pqLO.push(t);
                        HIcount++;
                    }
                    else{
                        t.discarded = true;
                        pendingJobs.push(t);
                        LOcount++;
                    }
                    // printf("Job %d of Task %d arrived with criticality %d\n", t.job_id + 1, t.task_id + 1, t.criticality);
                    cur_task_info[j]++;
                }
            }
            return HIcount;
            
        }

        
        void computeCurrentSlack2(double cur_time, int k, custom_priority_queue<JOB>& rq, JOB cur_job, double extra_time){
            custom_priority_queue<JOB> tmp;  // to store slackless-discarded jobs
            custom_priority_queue<JOB> rqCopy = rq;
            double disc_used_slack = 0.0;
            while(!discardedJobs.empty()){
                JOB disc_job = discardedJobs.top();
                discardedJobs.pop();
                if(disc_job.deadline <= cur_time){
                    dynamic_discarded++;
                    continue;
                }
                if(current_level[coreNumber] > k){
                    disc_job.deadline = disc_job.original_deadline; // reseting deadline to original deadline
                }
                double slack = disc_job.deadline - cur_time;
                // printf("\ndisc_job.deadline: %lf\n", disc_job.deadline);


                if(current_level[coreNumber] > k){
                    cur_job.deadline = cur_job.original_deadline; // resetting deadline to original deadline
                }
                // double exeFrac = min(1.0, (disc_job.deadline - cur_time)/(cur_job.deadline - cur_time));
                // slack -= exeFrac * (execution_time[cur_job.task_id][min(current_level[coreNumber], cur_job.criticality)-1] - cur_job.executed_time);
                // printf("%T%d(J%d):%lf, ", cur_job.task_id+1, cur_job.job_id+1, cur_job.deadline);

                while(!rqCopy.empty()){
                    JOB job = rqCopy.top();
                    rqCopy.pop();
                    double exeFrac = min(1.0, (disc_job.deadline - cur_time)/(job.deadline - cur_time)); // there would be some undiscarded job also in readyQueue
                    slack -= exeFrac * (execution_time[job.task_id][min(current_level[coreNumber], job.criticality)-1] - job.executed_time_hfq);
                    // printf("^T%d(J%d):%lf, ", job.task_id+1, job.job_id+1, job.deadline);
                }
                for(auto job:jobs){
                    if(job.arrival_time <= cur_time or job.arrival_time >= disc_job.deadline or job.criticality < current_level[coreNumber])
                        continue;
                    if(current_level[coreNumber] > k){
                        job.deadline = job.original_deadline; // resetting deadline to original deadline
                    }
                    double exeFrac = min(1.0, (disc_job.deadline - job.arrival_time)/(job.deadline - job.arrival_time));
                    slack -= exeFrac * execution_time[job.task_id][current_level[coreNumber]-1];
                    // printf("T%d(J%d):%lf, ", job.task_id+1, job.job_id+1, job.deadline);

                }
                printf("\n[Core %d] slack: %lf\n", coreNumber, slack);
                if((slack - disc_used_slack) >= (execution_time[disc_job.task_id][disc_job.criticality-1] - disc_job.executed_time_hfq)){
                    printf("\n[Core %d] T%d(J%d) with remain time:%lf undiscarded due to slack of %lf\n\n", coreNumber, disc_job.task_id + 1, disc_job.job_id + 1, execution_time[disc_job.task_id][disc_job.criticality-1] - disc_job.executed_time_hfq ,slack);
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
        void pushRqIntoDiscardedJobs(custom_priority_queue<JOB>& pqLO){ 
            custom_priority_queue<JOB> tmp;
            while(!pqLO.empty()){
                JOB cur_job = pqLO.top();
                pqLO.pop();
                if(cur_job.criticality < current_level[coreNumber]){
                    cur_job.discarded = true;
                    discardedJobs.push(cur_job);
                }
                else{
                    tmp.push(cur_job);
                }
            }
            pqLO = tmp;
        }

        void doSlackWorks(double& cur_time, int& k, custom_priority_queue<JOB>& pqLO, JOB cur_job, double extra_time){
            if(pendingJobs.empty())
                return;
            while(!pendingJobs.empty()){
                JOB job = pendingJobs.top();
                pendingJobs.pop();
                discardedJobs.push(job);
            }
            // computeCurrentSlack(cur_time, k, pqLO, extra_time);
            computeCurrentSlack2(cur_time, k, pqLO, cur_job, extra_time);
        }

        double calculateFrequency(){
            double sum = 0.0;
            for(int i = 0; i < num_of_task; i++){
                sum += utilTask[i];
            }
            // printf("sum = %lf\n", sum);
            int index = upper_bound(availableFrequency.begin(), availableFrequency.end(), sum) - availableFrequency.begin();
            return availableFrequency[min(index, (int)availableFrequency.size()-1)];
        }

        double computeCoreSlack(double cur_time, int k, int core, custom_priority_queue<JOB> rq, JOB cur_job, JOB avoid_job){
            double miniSlack = INT_MAX;
            custom_priority_queue<JOB> rqCopy = rq;
            if(cur_job.job_index != avoid_job.job_index)
                rqCopy.push(cur_job);
            vector<JOB> rqJobs;
            while(!rqCopy.empty()){
                JOB job = rqCopy.top();
                rqCopy.pop();
                rqJobs.push_back(job);
            }
            for(auto job:jobs){
                if(job.arrival_time <= cur_time or job.criticality < current_level[core] or job.job_index == avoid_job.job_index)
                    continue;
                if(current_level[core] > k){
                    job.deadline = job.original_deadline; // resetting deadline to original deadline
                }
                rqJobs.push_back(job);
            }
            for(int i = 0; i < rqJobs.size(); i++){ // QOS for already undiscarded jobs but not for those which will be discarded later and might have been undiscarded due to slack
                JOB job = rqJobs[i];
                if(current_level[core] > k){ // DEBUG
                    job.deadline = job.original_deadline; // resetting deadline to original deadline
                }
                double slack = job.deadline - cur_time;
                for(int j = 0; j < rqJobs.size(); j++){
                    JOB other_job = rqJobs[j]; // we need to consider execution time of current job also since we are finding slack of system due to this job and later checking if slack > 0
                    if(other_job.arrival_time >= job.deadline)
                        continue;
                    if(current_level[core] > k){ // DEBUG
                        other_job.deadline = other_job.original_deadline; // resetting deadline to original deadline
                    }
                    double exeFrac = min(1.0, (job.deadline - cur_time)/(other_job.deadline - cur_time)); // there would be some undiscarded job also in readyQueue
                    slack -= exeFrac * (execution_time[other_job.task_id][min(current_level[core], other_job.criticality)-1] - other_job.executed_time_hfq);
                    // printf("T%d(J%d):%lf, ", other_job.task_id+1, other_job.job_id+1, slack);
                }
                // 10 - 6 - 4
                miniSlack = min(miniSlack, slack);
                miniSlack = round(miniSlack * 1e6) / 1e6;
                // printf("miniSlack[T%d(J%d)]: %lf\n", job.task_id+1, job.job_id+1, miniSlack);
            }
            return miniSlack;
        }
        bool migratable(double cur_time, int migratedCore, JOB mig_job, int k){
            custom_priority_queue<JOB> coreRqCopy = rq[migratedCore];
            coreRqCopy.push(mig_job); // assuming migrating job inside ready queue of migrated core(eventhough it didn't arrive yet in real time)
            JOB dummy_job;
            dummy_job.job_index = -1;
            double slack = computeCoreSlack(cur_time, k, migratedCore, coreRqCopy, cur_jobs[migratedCore], dummy_job); // if u call doSlackWorks, then pendingJobsQueue empty is also considered
            double newFreq = calculateFrequencyOtherCore();
            if(slack > 0 && newFreq >= Foptimal_min && newFreq <= Foptimal_max)  // check for Foptimal of other core
                return true;
            return false;
        }

        void migrateJobs(double cur_time, int k){
            ////////////////////////////////////////////// migrate either upcoming job of every task or jobs in ready queue
            double maxSlack = INT_MIN;
            JOB jobToRemove;
            int typeJobToRemove = -1; // 0 - RQ job, 1 - upcoming job, 2 - current executing job
            custom_priority_queue<JOB> rqCopy = rq[coreNumber];
            vector<JOB> rqJobs;
            while(!rqCopy.empty()){
                JOB job = rqCopy.top();
                rqCopy.pop();
                rqJobs.push_back(job);
            }
            rqCopy = rq[coreNumber];
            for(auto job:rqJobs){
                rqCopy.remove(job);
                JOB dummy_job;
                dummy_job.job_index = -1;
                double slack = computeCoreSlack(cur_time, k, coreNumber, rqCopy, cur_jobs[coreNumber], dummy_job); // slack obtained by removing RQ job
                if (slack > maxSlack) {
                    maxSlack = slack;
                    jobToRemove = job;
                    typeJobToRemove = 0;
                }
                rqCopy.push(job); // Put back the previously removed job
            }

            int num_of_upcoming_checks = num_of_task;
            for(auto job:jobs){
                if(job.arrival_time <= cur_time)
                    continue;
                double slack = computeCoreSlack(cur_time, k, coreNumber, rqCopy, cur_jobs[coreNumber], job); // slack obtained by removing upcoming job
                if (slack > maxSlack) {
                    maxSlack = slack;
                    jobToRemove = job;
                    typeJobToRemove = 1;
                }
                num_of_upcoming_checks--;
                if(num_of_upcoming_checks == 0)
                    break;
            }
            double slack = computeCoreSlack(cur_time, k, coreNumber, rqCopy, cur_jobs[coreNumber], cur_jobs[coreNumber]); // slack obtained by removing current executing job in core
            if (slack > maxSlack) {
                maxSlack = slack;
                jobToRemove = cur_jobs[coreNumber];
                typeJobToRemove = 2;
            }

            ////////////////////////////////////////////////
            JOB mig_job;
            int task_id = jobToRemove.task_id;
            if(typeJobToRemove == 0){
                mig_job = jobToRemove;
                rq[coreNumber].remove(jobToRemove);
            }
            else if(typeJobToRemove == 1){
                mig_job = loadJob(task_id, (cur_task_info[task_id]) * period[task_id] + arrival_time[task_id]); // cur_time = (cur_task_info[i]) * period[i] + arrival_time[i]
            }
            else{
                mig_job = jobToRemove;
                cur_job_index = -1;
            }
            
            bool flag = false;
            for(int i = 0; i < NO_CORES; i++){
                if(i == coreNumber or shutdownable[coreNumber] or !shutdownable[i] or !migratable(cur_time, i, mig_job , k)) // check for LO criticality also
                    continue;
                migratedJobs[i].push_back(mig_job);
                flag = true;
                break;
            }
            if(flag){
                //
                //cur_task_info[earliest_deadline_task]++; // this job has arrived in migrated core at a later stage, hence update as arrived already in current core.
                migrateJobs(cur_time, k);
            }
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
            //custom_priority_queue<JOB> pqLO; // Min heap of pairs, where first is deadline and second is task id

            JOB cur_job;
            double cur_time = 0.0, prev_time = 0.0, extra_time = 0.0, freq = 1.0;
            

            double decision_point = findMinTask();
            if(decision_point != cur_time){
                total_idle_time += decision_point - cur_time;
                printf("[Core %d] %lf -> %lf => IDLE\n", coreNumber,cur_time, decision_point);
            }
            cur_time = decision_point;
            bool arrival = false;

            while(cur_time < time){
                bool idle = cur_job_index == -1;
                if(cur_job_index != -1){
                    cur_job.remain_time -= (cur_time - prev_time) * freq;
                    cur_job.executed_time_hfq += (cur_time - prev_time) * freq;
                    cur_job.executed_time += cur_time - prev_time;

                    cur_job.remain_time = round(cur_job.remain_time * 1e6) / 1e6;
                    cur_job.executed_time_hfq = round(cur_job.executed_time_hfq * 1e6) / 1e6;
                    cur_job.executed_time = round(cur_job.executed_time * 1e6) / 1e6;
                    current_jobs[cur_job.task_id] = cur_job;


                    if(abs(cur_job.remain_time) <= 0){  // remove completed job from CPU
                        // printf("Job %d of Task %d completed\n", cur_job.job_id + 1, cur_job.task_id + 1);                
                        utilTask[cur_job.task_id] = cur_job.executed_time_hfq / cur_job.deadline;
                        if(!cur_job.discarded){  


                            extra_time += execution_time[cur_job.task_id][current_level[coreNumber]-1] - cur_job.executed_time_hfq;
                            if(cur_job.executed_time < COMPLETION_LB * execution_time[cur_job.task_id][current_level[coreNumber]-1]){   // DOUBT
                                printf("[Core %d] T%d(J%d) completed early\n", coreNumber, cur_job.task_id + 1, cur_job.job_id + 1);
                                doSlackWorks(cur_time, k, rq[coreNumber], cur_job, extra_time);
                            }
                        }
                        cur_job_index = -1; 
                        prev_job_index = -1;
                    }

                    if(!cur_job.discarded && (cur_job.executed_time_hfq >= execution_time[cur_job.task_id][current_level[coreNumber]-1])){ // CHANGE >= to >
                        printf("\n[Core %d] Criticality Level changed from %d to %d\n\n", coreNumber, current_level[coreNumber], current_level[coreNumber]+1);
                        current_level[coreNumber]++;
                        printf("[Core %d]Job %d of Task %d preempted(%lf)\n", coreNumber, cur_job.job_id + 1, cur_job.task_id + 1, cur_job.remain_time);
                        cur_job.preempted_time = cur_time;  // current_level[coreNumber] changing is also a preemption point right now
                        rq[coreNumber].push(cur_job);
                        cur_job_index = -1;
                        pushRqIntoDiscardedJobs(rq[coreNumber]);
                        if(current_level[coreNumber] > k){
                            transferRq(rq[coreNumber]);
                        }
                        doSlackWorks(cur_time, k, rq[coreNumber], cur_job, extra_time);
                    }
                    
                    if(cur_job_index != -1){
                        // printf("Job %d of Task %d preempted(%lf)\n", cur_job.job_id + 1, cur_job.task_id + 1, cur_job.remain_time);
                        cur_job.preempted_time = cur_time;  // current_level[coreNumber] changing is also a preemption point right now
                        rq[coreNumber].push(cur_job);
                        cur_job_index = -1;
                    }
                }
                if(!arrival){ 
                    int HIcount = handleArrivals(cur_time, k, rq[coreNumber]);
                    if(HIcount > 0)
                        doSlackWorks(cur_time, k, rq[coreNumber], cur_job, extra_time);
                    arrival = true;  
                }
                
                
                
                if(!rq[coreNumber].empty()){
                    cur_job = rq[coreNumber].top();
                    rq[coreNumber].pop();

                    cur_job_index = cur_job.job_index;

                    if(cur_job.original_deadline <= cur_time){
                        printf("[Core %d] Job %d of Task %d missed deadline(%lf)\n", coreNumber, cur_job.job_id + 1, cur_job.task_id + 1, cur_job.deadline);
                        cur_job_index = -1;
                        continue;
                    }
                    if(cur_job.original_deadline - cur_time < cur_job.remain_time){
                        printf("[Core %d] Job %d of Task %d will miss deadline\n", coreNumber, cur_job.job_id + 1, cur_job.task_id + 1);
                        cur_job_index = -1;
                        continue;
                    }

                    if(prev_job_index != -1 && cur_job_index != prev_job_index){
                        preemptions++;
                        printf("[Core %d] Preemption at t = %lf\n", coreNumber, cur_time);
                    }
                    prev_job_index = cur_job_index;





                    utilTask[cur_job.task_id] = (1.0*execution_time[cur_job.task_id][criticality[cur_job.task_id]-1])/cur_job.deadline;
                    freq = calculateFrequency();
                    if(freq > Foptimal_max){
                        migrateJobs(cur_time , k); // new freq can go below Foptimal_min
                    }
                    freq = calculateFrequency();
                    if(freq < Foptimal_min){
                        acceptJobs(cur_time); // accept jobs from other cores UNTIL freq becomes within Foptimal
                    }
                    freq = calculateFrequency();

                    printf("\n[Core %d] Frequency changed to %lf\n", coreNumber, freq);
                    printf("[Core %d] Old remain time: %lf, ", coreNumber, cur_job.remain_time);
                    // cur_job.remain_time /= freq;
                    printf("[Core %d] New remain time: %lf\n", coreNumber, cur_job.remain_time / freq);

                    if(!cur_job.discarded)
                        decision_point = min({cur_time + cur_job.remain_time / freq, cur_time + execution_time[cur_job.task_id][current_level[coreNumber]-1] - cur_job.executed_time_hfq, findMinTask()});
                    else    
                        decision_point = min({cur_time + cur_job.remain_time / freq, findMinTask()});

                    
                    if(decision_point == cur_time + cur_job.remain_time){
                        if(!cur_job.discarded)
                            printf("[Core %d] %lf -> %lf => T%d(J%d)*\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                        else{
                            printf("[Core %d] %lf -> %lf => Discarded T%d(J%d)*\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                            // printf("cur_job.remain_time: %lf\n", cur_job.remain_time);
                            // printf("execution_time[cur_job.task_id][current_level[coreNumber]-1] = %lf\n", execution_time[cur_job.task_id][current_level[coreNumber]-1]);
                        }
                    }
                    else{
                        if(!cur_job.discarded)
                            printf("[Core %d] %lf -> %lf => T%d(J%d)\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                        else{
                            printf("[Core %d] %lf -> %lf => Discarded T%d(J%d)\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                            // printf("cur_job.remain_time: %lf\n", cur_job.remain_time);
                            // printf("execution_time[cur_job.task_id][current_level[coreNumber]-1] = %lf\n", execution_time[cur_job.task_id][current_level[coreNumber]-1]);
                        }
                    }
                }
                else{
                    decision_point = findMinTask();
                    total_idle_time += decision_point - cur_time;
                    printf("[Core %d] %lf -> %lf => IDLE\n", coreNumber, cur_time, decision_point);
                }
                usleep(SLEEP_TIME);
                prev_time = cur_time;
                cur_time = decision_point;
                arrival = false;

            }
        }

        void initUtilTask(){
            for(int i = 0; i < num_of_task; i++){
                utilTask[i] = (1.0*execution_time[i][criticality[i]-1])/deadline[i];
                printf("[Core %d] utilTask[%d] = %lf\n", coreNumber, i, utilTask[i]);
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
                printf("[Core %d] T%d(J%d) : Deadline = %lf\n", coreNumber, job.task_id+1, job.job_id+1, job.deadline);
            }
        }
        void edfvd(int hyperperiod) {
            int k = offlinePreprocessing(); //todo
            printf("k: %d\n", k);
            loadJobsInHyperperiod(hyperperiod);
            // printJobs();
            num_of_jobs = jobs.size();
            // slackTablePreCompute();
            // printSlackTable();
            initUtilTask();
            runtimeScheduling(k, hyperperiod);
        }
        void init(vector<int> myTasks)
        {
            for(int i = 0; i < num_of_task; i++)
            {
                global_id[i] = myTasks[i];
                arrival_time[i] = global_arrival_time[myTasks[i]];
                deadline[i] = global_deadline[myTasks[i]];
                criticality[i] = global_criticality[myTasks[i]];
                period[i] = global_period[myTasks[i]];
                virtualDeadlines[i] = global_virtualDeadlines[myTasks[i]];
                max_criticality = max(max_criticality, criticality[i]);
                for(int j = 0; j < criticality[i]; j++)
                {
                    execution_time[i][j] = global_execution_time[myTasks[i]][j];
                }
            }
            hyperperiod = findLCM();
        }
};
 
class ExceptionCore
{
    public:
        void print(int x)
        {
            cout << x << endl;
        }
};
void get_task_info(){
	FILE* file;
	file = fopen("tasks-vd.txt", "r");
	if(file == NULL){
		cout << "File not found" << endl;
		exit(0);
	}
	fscanf(file, "%d", &global_num_of_task);
	for(int i = 0; i < global_num_of_task; i++){
		fscanf(file, "%lf %d %lf %d", &global_arrival_time[i], &global_period[i], &global_deadline[i], &global_criticality[i]);
        for(int j = 0; j < global_criticality[i]; j++){
            fscanf(file, "%lf", &global_execution_time[i][j]);
        }
	}
}
int main()
{
    //create thread and pass 5 as parameter
    ShutDownableCore core, core1;
    vector<int> tasks;
    tasks.push_back(0);
    tasks.push_back(1);
    tasks.push_back(2);
    get_task_info();
    thread t2(&ShutDownableCore::mainEq, &core, tasks, 0);
    thread t1(&ShutDownableCore::mainEq, &core1, tasks, 1);
    t2.join();
    t1.join();
}