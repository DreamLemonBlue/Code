#include<iostream.h>
#include<math.h>
#include<stdlib.h>
#include<fstream.h>																						 //����ͷ


/*-------------------------------------------------------------------------------------*/
//�������޸���ľ���

double eps=1e-6;

/*-------------------------------------------------------------------------------------*/
//�������޸���ĺ���f(x)

double f(double x)
{
	return (4 * powl(x, 3) - 8 * powl(x, 2) + 2 * x + 1);
}

/*-------------------------------------------------------------------------------------*/
//�������޸���ĵ���df(x)

double df(double x)
{
	return (12 * powl(x, 2) - 16 * x + 2);
}

/*-------------------------------------------------------------------------------------*/
//�������޸���ĵ���ʽg(x)

double g(double x)
{
	double p=1.0/3;
	return (powl((8*powl(x,2)-2*x-1),p));
}

/*-------------------------------------------------------------------------------------*/
//���´��������޸�
//����������

void welcome();																								//��ӭҳ��

class solve;																										//������

void method_select(solve &a);																	//��������

double half(double a,double b,double eps);												//���ַ��㷨��������
double newn(double x0, double eps);															//ţ�ٵ������㷨��������
double itr(double (*fun)(double), double x0, double eps);						//�򵥵������㷨��������
double chord(double x0, double eps);														//�ҽط��㷨��������
double stef(double x0, double x);																//ʷ�ٷ����㷨�����������
double steffensen(double x0, double eps);												//ʷ�ٷ����㷨��������

double gc(double x0,double x)																		//�ҽط���ʽ
{
	return (x-f(x)/(f(x)-f(x0))*(x-x0));
}
double g1(double x)																						//�򵥵�����ʽ1
{
	return (-1.0) / (4 * pow(x, 2) - 8 * x + 2);
}
double g2(double x)																						//�򵥵�����ʽ2
{
	double t=1.0/2;
	return pow(((4*pow(x,3)+2*x+1)*(1.0/8)),t);
}
double g3(double x)																						//�򵥵�����ʽ3,������ĵ���ʽg(x)
{
	double t=1.0/3;
	return pow(((8*pow(x,2)-2*x-1)*(1.0/4)),t);
}
double g4(double x0)																					//ţ�ٵ������ƽ�
{
	return (x0 - f(x0) / df(x0));
}
double g5(double x)																						//ʷ�ٷ����㷨������ʽ
{
	return ((-4 * pow(x, 3) + 8 * pow(x, 2) - 1) * (1.0 / 2));
}

class solve																										//�ⷨ��
{
	double a,b;																									//��������˵�
	double h;																									//����ֵ
	double c[5][2];																							//����˵㴢��
	int i;																												//��������
	double result[3];																						//�������
public:
	solve(double p=0,double q=1,double r=0.1)												//���캯��
	{
		a=p;
		b=q;
		h=r;
	}
	void filter()																									//ɸѡ���䲢����
	{
		i=0;
		double a0,b0;
		a0=a;
		b0=a+h;
		do
		{
			if(f(a0)*f(b0)<0)
			{
				c[i][0]=a0;
				c[i][1]=b0;
				i++;
			}
			a0=b0;
			b0+=h;	
		}while(b0<=b);
	}
	void halfs()																									//���ַ�
	{
		filter();
		double x;
		cout<<"IT\ta1\tb1\tx\tf(x)\n";
		for(int j=0;j<i;j++)
		{
			a=c[j][0];
			b=c[j][1];
			x=half(a,b,eps);
			result[j] = x;
			cout<<"x="<<x<<endl;
		}
	}
	void itrs()																									//�򵥵�����
	{
		filter();
		double x;
		cout << "IT\ta1\tb1\tx\tf(x)\n";
		x = itr(g1, c[0][0], eps);
		result[0] = x;
		cout << "x1=" << x << endl;
		x = itr(g2, c[1][0], eps);
		result[1] = x;
		cout << "x2=" << x << endl;
		x = itr(g3, c[2][0], eps);
		cout << "x3=" << x << endl;
		result[2]=x;
	}
	void newns()																								//ţ�ٵ�����
	{
		filter();
		double x;
		for (int j = 0; j < i; j++)
		{
			cout << "IT\tx0\t\tx\n";
			x = newn(c[j][0], eps);
			result[j] = x;
			cout << "\nx=" << x << endl;
		}
	}
	void chords()																								//�ҽط�
	{
		filter();
		double x;
		cout << "IT\tx0\t\tx\n";
		for (int q = 0; q < i; q++)
		{
			x = chord(c[q][0], eps);
			result[q] = x;
			cout << "\nx=" << x << endl;
		}
	}
	void steffensens()																						//Steffensen��
	{
		filter();
		double x;
		cout << "IT\tx0\t\tx\n";
		for (int r = 0; r < i; r++)
		{
			x = steffensen(c[r][0], eps);
			result[r] = x;
			cout << "x=" << x << endl;
		}
	}
	void output()																								//������ļ�
	{
		ofstream outfile("answer.txt");
		for(int i=0;i<3;i++)
		{
			outfile<<"x"<<i<<"="<<result[i]<<'\n';
		}
		outfile.close();
	}
	
};

void main()																										//������
{
	welcome();
	double m,n;																								//����˵�
	cout<<"����������[a,b]�˵�:"<<endl;
	cout<<"a=";
	cin>>m;
	cout<<"b=";
	cin>>n;
	solve a(m,n);
	method_select(a);
	a.output();
}

void method_select(solve &a)
{
	system("cls");
	cout<<"����������ѡ�񷽷���Ӧ������:"<<endl;
	cout<<"1:���ַ�\n2:�򵥵�����\n3:Steffensen��\n4:ţ�ٷ�\n5:�ҽط�\n"<<endl;
	int n;
	cin>>n;
	switch(n)//��֧
	{
	case 1:		a.halfs();
		break;
	case 2:		a.itrs();
		break;
	case 3:		a.steffensens();
		break;
	case 4:		a.newns();
		break;
	case 5:		a.chords();
		break;
	}
}

void welcome()
{
	cout<<"��ӭ��ʹ�ñ�����!\n";
	cout<<"�������������ţ��ں��֣���������ͬ���\n";
	cout<<"�޸ķ�������Դ�ļ����޸ĺ������������͵�������\n";
	cout<<"Ĭ�Ϸ���:4x^3-8x^2+2x+1=0\n";
	cout<<"ѡ�񷽷�������������,�õ��Ľ�ᱻ������answer.txt��\n";
}
//�����ǲ�ͬ�㷨���Եĺ���
double half(double a,double b,double eps)
{
	int it=0;
	double x,a1,b1,y0,y;
	a1=a;
	b1=b;
	y0=f(a);
	do
	{
		x=(a1+b1)/2;
		cout<<it<<'\t'<<a1<<'\t'<<b1<<'\t'<<x<<'\t'<<f(x)<<endl;
		y=f(x);
		if(y*y0>0)
		{
			a1=x;
			y0=f(a1);
		}
		else
		{
			b1=x;		
		}
		it++;
	}
	while(fabs(y)>eps);
	return(x);
}//���ַ�
double newn(double x0, double eps)//ţ�ٵ�����
{
	int it = 0, e;
	double x;
	do
	{
		x = g4(x0);
		cout << it << '\t' << x0 << '\t' << '\t' << x << endl;
		e = fabs((x - x0) / x) > eps;
		if (e) { x0 = x; it++; }
	} while (e);
	return(x);
}
double itr(double (*fun)(double), double x0, double eps)//�򵥵�����
{
	int it = 0, e;
	double x;
	do
	{
		x = fun(x0);
		cout << it << '\t' << x0 << '\t' << '\t' << x << endl;
		e = fabs((x - x0) / x) > eps;
		if (e) { x0 = x; it++; }
	} while (e);
	return(x);
}
double chord(double x0, double eps)//�ҽط�
{
	int it = 0, e;
	double x, x1;
	if (f(x0) > 0) x1 = x0 + 0.1;
	else x1 = x0 - 0.1;
	do
	{
		x = gc(x0, x1);
		cout << it << '\t' << x0 << '\t' << '\t' << x << endl;
		e = fabs((x1 - x0) / x1) > eps;
		if (e) { x0 = x1; x1 = x; it++; }
	} while (e);
	return(x);
}
double stef(double x0, double x)//Steffensen�㷨���
{
	double y, z;
	y = g5(x0);
	z = g5(y);
	return (x0 - (powl((y - x0), 2)) / (z - 2 * y + x0));
}
double steffensen(double x0, double eps)//Stenffensen��
{
	int it = 0, e;
	double x, x1;
	if (f(x0) > 0)
		x1 = x0 + 0.1;
	else
		x1 = x0 - 0.1;
	do
	{
		x = stef(x0, x1);
		cout << it << '\t' << x0 << '\t' << '\t' << x << endl;
		e = fabs((x1 - x0) / x1) > eps;
		if (e)
		{
			x0 = x1;
			x1 = x;
			it++;
		}
	} while (e);
	return(x);
}