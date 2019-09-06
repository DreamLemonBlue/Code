/*
Author:DreamLemonBlue
程序说明:解线性方程组(唯一解)
*/

#include<iostream>
#include<cmath>
#include<fstream>
#include<cstdlib>//头文件
using namespace std;
#define N 20//方程阶数上限
#define M 500//迭代次数上限，可修改，超过视为迭代发散

class leg;//声明类
void welcome();//欢迎页面与初始页面
int input_select(int &num, double xn[], double ar[][N]);//选择输入方法(控制台/文件流)分支以及输入
void method_select(leg &a);//选择解线性方程组的方法

class leg//linear equation group,线性方程组
{
	int n;//方程组中方程个数，或系数矩阵阶数
	double martix[N][N];//记录系数矩阵
	double p[N];//记录常数项
	double x[N];//解的储存
	double d;//系数行列式
public:
	leg(double a[][N], double b[], int num)//构造函数
	{
		n = num;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				martix[i][j] = a[i][j];
			}
		for (int k = 0; k < n; k++)
		{
			p[k] = b[k];
		}
		d = determinant(martix,n);
	}
	double determinant(double arcs[][N],int n0)//求行列式
	{
		if (n0 == 1)
		{
			return arcs[0][0];
		}
		double ans=0.0;
		double temp[N][N] = { 0.0 };
		int i, j, k;
		for (i = 0; i < n0; i++)
		{
			for (j = 0; j < n0 - 1; j++)
			{
				for (k = 0; k < n0 - 1; k++)
				{
					temp[j][k] = arcs[j + 1][(k >= i) ? k + 1 : k];
				}
			}
			double t = determinant(temp, n0 - 1);
			if (i % 2 == 0)
			{
				ans += arcs[0][i] * t;
			}
			else
			{
				ans -= arcs[0][i] * t;
			}
		}
		return ans;
	}
	bool check()//检查线性方程组是否有唯一解
	{
		if (!d)
		{
			cout << "系数方阵不满秩" << endl;
			cout << "输入的方程组没有唯一解" << endl;//不满秩退出
			return false;
		}
		else
			return true;
	}
	void Cramer()//克拉默法求解
	{
		double dm[N];
		double temp[N][N];
		for (int i = 0; i < n; i++)
		{
			for (int x = 0; x < n; x++)
			{
				for (int y = 0; y < n; y++)
				{
					temp[x][y] = martix[x][y];
				}
			}//初始化为martix
			for (int z = 0; z < n; z++)
			{
				temp[z][i] = p[z];
			}
			dm[i] = determinant(temp, n);//替代一列
		}
		for (int m = 0; m < n; m++)
		{
			x[m] = dm[m] / d;
		}
		output();
	}
	void Gauss()//高斯列主元消去法求解
	{
		double temp[N][N];
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				temp[i][j] = martix[i][j];
			}
		}
		for (int r = 0; r < n; r++)
		{
			temp[r][n] = p[r];
		}//将类中数据改成增广矩阵形式
		ColPivot(temp,n,x);
		output();
	}
	void Gauss_Seidel()//高斯-赛德尔迭代法求解
	{
		double temp[N][N];
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				temp[i][j] = martix[i][j];
			}
		}
		for (int r = 0; r < n; r++)
		{
			temp[r][n] = p[r];
		}//同理，将类中数据改成增广矩阵形式
		GS(temp,n,x);
	}
	void output()//分别在控制台和文件流中输出线性方程组的解
	{
		cout << "此线性方程组的解如下:\n";
		for (int i = 0; i < n; i++)
		{
			cout << "x[" << i << "]=" << x[i] << endl;
		}
		cout << "答案已保存至\"answer.txt\"\n";
		ofstream outfile("answer.txt");
		outfile << "输入的方程:" << endl;
		for (int t = 0; t < n; t++)
		{
			outfile << martix[t][0] << "x1";
			for (int y = 1; y < n; y++)
			{
				if (martix[t][y] > 0)
					outfile << "+" << martix[t][y] << "x" << y + 1;
				if (martix[t][y] < 0)
					outfile << martix[t][y] << "x" << y + 1;
				if (y == (n - 1))
					outfile << "=" << p[t] << endl;
			}
		}
		outfile << "\n\n";
		outfile << "此线性方程组的解如下:\n";
		for (int j = 0; j < n; j++)
		{
			outfile<< "x[" << j << "]=" << x[j] << endl;
		}

	}//输出到控制台以及文件
	void showinf()//显示输入的线性方程组的信息
	{
		cout << "输入的方程:" << endl;
		for (int x = 0; x < n; x++)
		{
			cout << martix[x][0] << "x1";
			for (int y = 1; y < n; y++)
			{
				if (martix[x][y] > 0)
					cout << "+" << martix[x][y] << "x" << y + 1;
				if (martix[x][y] < 0)
					cout << martix[x][y] << "x" << y + 1;
				if (y == (n - 1))
					cout << "=" << p[x] << endl;
			}
		}
		cout << "\n\n";
	}
	void ColPivot(double c[][N], int n, double x[])//高斯主元素消去算法
	{
		int i, j, t, k;
		double p;
		for (i = 0; i < n - 1; i++)
		{
			k = i;
			for (j = i + 1; j < n; j++)
			{
				if ((fabs(c[j][i])) < (fabs(c[k][i])))
					k = j;
				if (k != i)
				{
					for (j = i; j < n + 1; j++)
					{
						p = c[i][j];
						c[i][j] = c[k][j];
						c[k][j] = p;
					}
				}
			}
			for (j = i + 1; j < n; j++)
			{
				p = (c[j][i]) / (c[i][i]);
				for (t = i; t < n + 1; t++)
				{
					c[j][t] -= p * (c[i][t]);
				}
			}
		}
		for (i = n - 1; i > -1; i--)
		{
			for (j = n - 1; j > i; j--)
			{
				(c[i][n]) -= x[j] * (c[i][j]);
			}
			x[i] = (c[i][n]) / (c[i][i]);
		}
	}
	void GS(double a[][N], int n, double x[])//Gauss-Seidel迭代算法
	{
		int i, j, k = 1;
		double d, dx, eps;
		for (i = 0; i <= n - 1; i++)
		{
			x[i] = 0.0;
		}
		eps = 0;
		while (1)
		{
			for (i = 0; i <= n - 1; i++)
			{
				d = 0;
				for (j = 0; j <= n - 1; j++)
				{
					if (j == i)
						continue;
					d += a[i][j] * x[j];
				}
				dx = ((a[i][n] - d) / (a[i][i]));
				eps = fabs(dx - x[i]);
				x[i] = dx;
			}
			if (eps < 1e-6)
			{
				cout << "最后一次迭代差值为:" << eps << endl;
				cout << "迭代收敛!迭代次数为:" << k << endl;
				output();
				return;
			}
			if (k > M)
			{
				cout << "最后一次迭代差值为:" << eps << endl;
				cout << "迭代发散!迭代次数为:" << M << endl;
				return;
			}
			k++;
		}
	}
};

int main()//主函数
{
	int loop=1;
	do
	{
		int n;
		double ar[N][N];
		double xn[N];
		welcome();
		if (input_select(n, xn, ar))
		{
			leg a(ar, xn, n);
			a.showinf();
			if (a.check())
			{
				method_select(a);
			}
		}
		cout << "\n输入任意非0值以返回程序初始页面或输入0结束本程序\n";
		cin >> loop;
		system("cls");
	} while (loop);//循环实现返回初始菜单
	return 0;
}

void welcome()
{
	cout << "欢迎你使用本程序\n"
			<< "本程序用于解线性方程组,如果有唯一解会将结果输出\n"
			<< "默认设置最大阶数为20\n";
}
int input_select(int &num,double xn[],double ar[][N])
{
	int method;
	cout << "\n\n1:在控制台输入线性方程组\n"
		<< "2:在文件流中输入线性方程组\n"
		<< "\n输入对应的数字以选择输入方法:";
	cin >> method;
	system("cls");
	switch (method)
	{
	case 1: {cout << "你选择了控制台输入\n\n";
		cout << "请输入方程阶数n=";
		cin >> num;
		cout << "请输入系数方阵\n";
		for (int i = 0; i < num; i++)
		{
			for (int j = 0; j < num; j++)
			{
				cin >> ar[i][j];
			}
		}
		cout << "\n系数方阵输入完毕\n";
		cout << "\n请输入常数项\n";
		for (int k = 0; k < num; k++)
		{
			cin >> xn[k];
		}
		cout << "\n输入完毕\n";
		return 1; }
	case 2:{cout << "你选择了文件流输入\n\n";
		cout << "请按以下格式输入:\n";
		cout << "[方程阶数]\n";
		cout << "[空行]\n";
		cout << "[系数方阵]\n";
		cout << "[空行]\n";
		cout << "[常数项]\n\n";
		cout << "请将文件保存为\"input.txt\"";
		int flag = 0;
		cout << "输入完成请输入1,中止输入请输入0\n";
		cin >> flag;
		if (flag)
		{
			cout << "确认输入完成\n";
			ifstream infile("input.txt");
			infile >> num;
			for (int i = 0; i < num; i++)
			{
				for (int j = 0; j < num; j++)
				{
					infile >> ar[i][j];
				}
			}
			for (int k = 0; k < num; k++)
			{
				infile >> xn[k];
			}
			infile.close();
			return 1;
		}
		else
		{
			cout << "输入中止!\n";
			return 0;
		}
	}
	default:		cout << "你输入的数字不正确\n";
						return 0;
	}

}
void method_select(leg &a)
{
	system("cls");
	int m;
	a.showinf();
	cout << "1:Cramer法\n"<< "2:高斯主元素消去法\n"<< "3:Gauss-Seidel迭代法\n"<< "请输入选择的方法对应的数字:";
	cin >> m;
	cout << endl;
	switch (m)
	{
	case 1:			cout << "你选择了Cramer法\n";
						a.Cramer();
						break;
	case 2:			cout << "你选择了高斯主元素消去法\n";
						a.Gauss();
						break;
	case 3:			cout << "你选择了Gauss-Seidel迭代法\n";
						a.Gauss_Seidel();
						break;
	}
}