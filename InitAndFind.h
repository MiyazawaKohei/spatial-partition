using namespace std;
vector<double*> Init_Vectors_inB(int n_b);
void Close_Vectors(vector<double*> B);
vector<double> Create_Random_Normarized_Vector();
double L2_Norm(vector<double> v);
vector<double> Normarize_Vector(vector<double> v);
vector<vector<double>> Create_Random_Rotation_Matrix();
int MyRandom_Int(int i);
vector<int> Create_Random_Permutaiton(int n_b);
vector<double*> Init_Vectors_inA(vector<double*> B,vector<vector<double>> rotation, vector<int> permutation,int n_a, double SDofnoise);
void Fix_Gravity_Point_to_Origin(vector<double*> vectors);
Eigen::Matrix<double,3,3> Compute_BAt(const vector<double*> A,const vector<double*> B);
double Maximize_trRBAt(const vector<double*> A,const vector<double*> B);
double Minimize_RMSD(const vector<double*> A,const vector<double*> B, const vector<int> sigma_i);
void make_initial_matrix(const vector<double*> A,const vector<double*> B,vector<vector<double>> *inequalities);

