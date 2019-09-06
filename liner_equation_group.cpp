/*
Author:DreamLemonBlue
����˵��:�����Է�����(Ψһ��)
*/

#include<iostream>
#include<cmath>
#include<fstream>
#include<cstdlib>//ͷ�ļ�
using namespace std;
#define N 20//���̽�������
#define M 500//�����������ޣ����޸ģ�������Ϊ������ɢ

class leg;//������
void welcome();//��ӭҳ�����ʼҳ��
int input_select(int &num, double xn[], double ar[][N]);//ѡ�����뷽��(����̨/�ļ���)��֧�Լ�����
void method_select(leg &a);//ѡ������Է�����ķ���

class leg//linear equation group,���Է�����
{
	int n;//�������з��̸�������ϵ���������
	double martix[N][N];//��¼ϵ������
	double p[N];//��¼������
	double x[N];//��Ĵ���
	double d;//ϵ������ʽ
public:
	leg(double a[][N], double b[], int num)//���캯��
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
	double determinant(double arcs[][N],int n0)//������ʽ
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
	bool check()//������Է������Ƿ���Ψһ��
	{
		if (!d)
		{
			cout << "ϵ����������" << endl;
			cout << "����ķ�����û��Ψһ��" << endl;//�������˳�
			return false;
		}
		else
			return true;
	}
	void Cramer()//����Ĭ�����
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
			}//��ʼ��Ϊmartix
			for (int z = 0; z < n; z++)
			{
				temp[z][i] = p[z];
			}
			dm[i] = determinant(temp, n);//���һ��
		}
		for (int m = 0; m < n; m++)
		{
			x[m] = dm[m] / d;
		}
		output();
	}
	void Gauss()//��˹����Ԫ��ȥ�����
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
		}//���������ݸĳ����������ʽ
		ColPivot(temp,n,x);
		output();
	}
	void Gauss_Seidel()//��˹-���¶����������
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
		}//ͬ�����������ݸĳ����������ʽ
		GS(temp,n,x);
	}
	void output()//�ֱ��ڿ���̨���ļ�����������Է�����Ľ�
	{
		cout << "�����Է�����Ľ�����:\n";
		for (int i = 0; i < n; i++)
		{
			cout << "x[" << i << "]=" << x[i] << endl;
		}
		cout << "���ѱ�����\"answer.txt\"\n";
		ofstream outfile("answer.txt");
		outfile << "����ķ���:" << endl;
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
		outfile << "�����Է�����Ľ�����:\n";
		for (int j = 0; j < n; j++)
		{
			outfile<< "x[" << j << "]=" << x[j] << endl;
		}

	}//���������̨�Լ��ļ�
	void showinf()//��ʾ��������Է��������Ϣ
	{
		cout << "����ķ���:" << endl;
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
	void ColPivot(double c[][N], int n, double x[])//��˹��Ԫ����ȥ�㷨
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
	void GS(double a[][N], int n, double x[])//Gauss-Seidel�����㷨
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
				cout << "���һ�ε�����ֵΪ:" << eps << endl;
				cout << "��������!��������Ϊ:" << k << endl;
				output();
				return;
			}
			if (k > M)
			{
				cout << "���һ�ε�����ֵΪ:" << eps << endl;
				cout << "������ɢ!��������Ϊ:" << M << endl;
				return;
			}
			k++;
		}
	}
};

int main()//������
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
		cout << "\n���������0ֵ�Է��س����ʼҳ�������0����������\n";
		cin >> loop;
		system("cls");
	} while (loop);//ѭ��ʵ�ַ��س�ʼ�˵�
	return 0;
}

void welcome()
{
	cout << "��ӭ��ʹ�ñ�����\n"
			<< "���������ڽ����Է�����,�����Ψһ��Ὣ������\n"
			<< "Ĭ������������Ϊ20\n";
}
int input_select(int &num,double xn[],double ar[][N])
{
	int method;
	cout << "\n\n1:�ڿ���̨�������Է�����\n"
		<< "2:���ļ������������Է�����\n"
		<< "\n�����Ӧ��������ѡ�����뷽��:";
	cin >> method;
	system("cls");
	switch (method)
	{
	case 1: {cout << "��ѡ���˿���̨����\n\n";
		cout << "�����뷽�̽���n=";
		cin >> num;
		cout << "������ϵ������\n";
		for (int i = 0; i < num; i++)
		{
			for (int j = 0; j < num; j++)
			{
				cin >> ar[i][j];
			}
		}
		cout << "\nϵ�������������\n";
		cout << "\n�����볣����\n";
		for (int k = 0; k < num; k++)
		{
			cin >> xn[k];
		}
		cout << "\n�������\n";
		return 1; }
	case 2:{cout << "��ѡ�����ļ�������\n\n";
		cout << "�밴���¸�ʽ����:\n";
		cout << "[���̽���]\n";
		cout << "[����]\n";
		cout << "[ϵ������]\n";
		cout << "[����]\n";
		cout << "[������]\n\n";
		cout << "�뽫�ļ�����Ϊ\"input.txt\"";
		int flag = 0;
		cout << "�������������1,��ֹ����������0\n";
		cin >> flag;
		if (flag)
		{
			cout << "ȷ���������\n";
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
			cout << "������ֹ!\n";
			return 0;
		}
	}
	default:		cout << "����������ֲ���ȷ\n";
						return 0;
	}

}
void method_select(leg &a)
{
	system("cls");
	int m;
	a.showinf();
	cout << "1:Cramer��\n"<< "2:��˹��Ԫ����ȥ��\n"<< "3:Gauss-Seidel������\n"<< "������ѡ��ķ�����Ӧ������:";
	cin >> m;
	cout << endl;
	switch (m)
	{
	case 1:			cout << "��ѡ����Cramer��\n";
						a.Cramer();
						break;
	case 2:			cout << "��ѡ���˸�˹��Ԫ����ȥ��\n";
						a.Gauss();
						break;
	case 3:			cout << "��ѡ����Gauss-Seidel������\n";
						a.Gauss_Seidel();
						break;
	}
}