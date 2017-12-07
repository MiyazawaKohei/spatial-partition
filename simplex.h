using namespace std;

void show_matrix(vector<vector<double>> *A,vector<int> *base);
void show_matrix(vector<vector<double>> *A);
int pivot(vector<vector<double>> *A, vector<int> *base);
int simplexinenr(vector<vector<double>> *A, vector<int> *base);
int simplex(vector<vector<double>> _A);