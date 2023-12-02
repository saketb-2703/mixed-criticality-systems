#include <bits/stdc++.h>
#include<unistd.h>  
using namespace std;

#define MAX_TASK 10
#define MAX_JOBS 100
#define NO_CORES 5
#define MAX_CRITICALITY 10
#define SLEEP_TIME 10000
#define COMPLETION_LB 0.8
#define PROCRAST_THRESHOLD 0.5
struct job{
	double arrival_time;
	int period;
	double executed_time;
    double executed_time_hfq;
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
int global_num_of_task;
double global_arrival_time[MAX_TASK], global_execution_time[MAX_TASK][MAX_CRITICALITY], global_deadline[MAX_TASK], global_virtualDeadlines[MAX_TASK], global_slackTable[MAX_CRITICALITY][MAX_JOBS][MAX_JOBS];
int global_criticality[MAX_TASK], global_period[MAX_TASK];

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
        int current_level = 1;
        int dynamic_discarded = 0;
        typedef struct job JOB;
        vector<double> utilTask{vector<double>(MAX_TASK,0.0)};
        vector<double> availableFrequency = {0.5, 0.75, 1.0};
        vector<JOB> jobs;
        double total_idle_time = 0.0, total_stolen_time = 0.0, slackSystem = 0.0;
        vector<JOB> current_jobs{vector<JOB>(MAX_TASK)};
        priority_queue<JOB> discardedJobs, pendingJobs;
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

        int handleArrivals(double cur_time, int k, priority_queue<JOB>& pqLO){
            int HIcount = 0, LOcount = 0;
            for (int j = 0; j < num_of_task; j++) {
                if (fmod((cur_time - arrival_time[j]), period[j]) == 0) {
                    JOB t = loadJob(j, cur_time);
                    if(t.criticality >= current_level){
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
                    slack -= exeFrac * (execution_time[other_job.task_id][min(current_level, other_job.criticality)-1] - other_job.executed_time_hfq);
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
                // printf("\ndisc_job.deadline: %lf\n", disc_job.deadline);


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
                    slack -= exeFrac * (execution_time[job.task_id][min(current_level, job.criticality)-1] - job.executed_time_hfq);
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
            double cur_time = 0.0, prev_time = 0.0, extra_time = 0.0, freq = 1.0;
            int cur_job_index = -1, prev_job_index = -1; // no job in CPU

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


                            extra_time += execution_time[cur_job.task_id][current_level-1] - cur_job.executed_time_hfq;
                            if(cur_job.executed_time < COMPLETION_LB * execution_time[cur_job.task_id][current_level-1]){   // DOUBT
                                printf("[Core %d] T%d(J%d) completed early\n", coreNumber, cur_job.task_id + 1, cur_job.job_id + 1);
                                doSlackWorks(cur_time, k, pqLO, cur_job, extra_time);
                            }
                        }
                        if(pqLO.empty()){
                            migrateJobs();
                            double procrast_slack = dpm(cur_time, pqLO, k);
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

                    if(!cur_job.discarded && (cur_job.executed_time_hfq >= execution_time[cur_job.task_id][current_level-1])){ // CHANGE >= to >
                        printf("\n[Core %d] Criticality Level changed from %d to %d\n\n", coreNumber, current_level, current_level+1);
                        current_level++;
                        printf("[Core %d]Job %d of Task %d preempted(%lf)\n", coreNumber, cur_job.job_id + 1, cur_job.task_id + 1, cur_job.remain_time);
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
                    int HIcount = handleArrivals(cur_time, k, pqLO);
                    if(HIcount > 0)
                        doSlackWorks(cur_time, k, pqLO, cur_job, extra_time);
                    if(idle){
                        double procrast_slack = dpm(cur_time, pqLO, k);
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
                
                
                
                if(!pqLO.empty()){
                    cur_job = pqLO.top();
                    pqLO.pop();

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
                        decision_point = min({cur_time + cur_job.remain_time / freq, cur_time + execution_time[cur_job.task_id][current_level-1] - cur_job.executed_time_hfq, findMinTask()});
                    else    
                        decision_point = min({cur_time + cur_job.remain_time / freq, findMinTask()});

                    
                    if(decision_point == cur_time + cur_job.remain_time){
                        if(!cur_job.discarded)
                            printf("[Core %d] %lf -> %lf => T%d(J%d)*\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                        else{
                            printf("[Core %d] %lf -> %lf => Discarded T%d(J%d)*\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                            // printf("cur_job.remain_time: %lf\n", cur_job.remain_time);
                            // printf("execution_time[cur_job.task_id][current_level-1] = %lf\n", execution_time[cur_job.task_id][current_level-1]);
                        }
                    }
                    else{
                        if(!cur_job.discarded)
                            printf("[Core %d] %lf -> %lf => T%d(J%d)\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                        else{
                            printf("[Core %d] %lf -> %lf => Discarded T%d(J%d)\n", coreNumber, cur_time, decision_point, cur_job.task_id + 1, cur_job.job_id + 1);
                            // printf("cur_job.remain_time: %lf\n", cur_job.remain_time);
                            // printf("execution_time[cur_job.task_id][current_level-1] = %lf\n", execution_time[cur_job.task_id][current_level-1]);
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
            int k = offlinePreprocessing();
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

class NonShutDownableCore
{
    public:
        void print(int x)
        {
            cout << x << endl;
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