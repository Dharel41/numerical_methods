#include<iostream>
#include<conio.h>
#include<fstream>
#include<ctime>
using namespace std;
fstream plik;
fstream plik2;
fstream plik3;
fstream results;
struct punkt
{
	double x, y;
};
double lagrangea()
{
	double a, b, suma = 0;

	cout << "Podaj ilosc znanych punktow" << endl;
	cin >> a;
	punkt *tab = new punkt[a];
	plik.open("plik.txt"), ios::in|ios::out;
	{
		cout << "Znane punkty :" << endl;
		if (plik.good())
		{
			int txtx, txty;
			for (int i = 0; i < a; i++)
			{
				plik >> txtx >> txty;
				tab[i].x = txtx;
				tab[i].y = txty;
				cout << "x:"<<txtx << " " << "y:" << txty << endl;
			}
		}
	}
	plik.close();

	double *tabl = new double[a];

	cout << "W ktorym punkcie chcesz znac wartosc?" << endl;
	cin >> b;

	for (int i = 0; i < a; i++)
	{
		tabl[i] = 1;
		for (int j = 0; j < a; j++)
		{
			if(i==j)
			{
				continue;
			}
			tabl[i] = tabl[i] * (b - tab[j].x) / (tab[i].x - tab[j].x);
		}
	}
	for (int i = 0; i < a; i++)
	{
		suma = suma + tab[i].y*tabl[i];
	}
	cout << "Wartosc w punkcie " << b << " wynosi " << suma;
	return suma;
}

double newton()
{
	int a, b;
	double suma, iloczyn;
	cout << "Podaj ilosc znanych punktow" << endl;
	cin >> a;
	punkt *tab = new punkt[a];
	cout << "Znane punkty :" << endl;
	plik.open("plik.txt"), ios::out | ios::in;
	{
		int txtx, txty;
		for (int i = 0; i < a; i++)
		{
			plik >> txtx >> txty;
			tab[i].x = txtx;
			tab[i].y = txty;
			cout << "x:" << txtx << " " << "y:" << txty << endl;
		}
	}
	plik.close();
	double **rzedy = new double*[a - 1];

	for (int i = 0; i < a - 1; i++)
	{
		rzedy[i] = new double[a - 1 - i];
	}

	for (int i = 0; i < a - 1; i++)
	{
		rzedy[0][i] = (tab[i + 1].y - tab[i].y) / (tab[i + 1].x - tab[i].x);
	}
	
	for (int j = 0; j < a - 2; j++)
	{
		for (int i = 0; i < a - 2 - j; i++)
		{
			rzedy[1 + j][i] = (rzedy[j][i + 1] - rzedy[j][i]) / (tab[i + 2 + j].x - tab[i].x);
		}
	}
	cout << "W ktorym punkcie chcesz znac wartosc?" << endl;
	cin >> b;
	suma = tab[0].y;
	for (int i = 0; i < a - 1; i++)
	{
		iloczyn = rzedy[i][0];
		for (int j = 0; j < i+1; j++)
		{
			iloczyn = iloczyn*(b-tab[j].x);
		}
		suma = suma + iloczyn;
	}
	cout << "Wartosc w punkcie " << b << " wynosi " << suma;
	return suma;
}

void wielomianowa2()
{
	double tab[3][4];
	plik.open("plik.txt"), ios::out | ios::in;
	{
		double txtx, w = 0, wx, wy, wz;
		for (int j = 0; j < 3; j++)
		{
			for (int i = 0; i < 4; i++)
			{
				plik >> txtx;
				tab[j][i] = txtx;
			}
		}
		w = tab[0][0] * tab[1][1] * tab[2][2] + tab[0][1] * tab[1][2] * tab[2][0] + tab[0][2] * tab[1][0] * tab[2][1] - tab[0][1] * tab[1][0] * tab[2][2] - tab[0][0] * tab[1][2] * tab[2][1] - tab[0][2] * tab[1][1] * tab[2][0];
		wx = tab[0][3] * tab[1][1] * tab[2][2] + tab[0][1] * tab[1][2] * tab[2][3] + tab[0][2] * tab[1][3] * tab[2][1] - tab[0][1] * tab[1][3] * tab[2][2] - tab[0][3] * tab[1][2] * tab[2][1] - tab[0][2] * tab[1][1] * tab[2][3];
		wy = tab[0][0] * tab[1][3] * tab[2][2] + tab[0][3] * tab[1][2] * tab[2][0] + tab[0][2] * tab[1][0] * tab[2][3] - tab[0][3] * tab[1][0] * tab[2][2] - tab[0][0] * tab[1][2] * tab[2][3] - tab[0][2] * tab[1][3] * tab[2][0];
		wz = tab[0][0] * tab[1][1] * tab[2][3] + tab[0][1] * tab[1][3] * tab[2][0] + tab[0][3] * tab[1][0] * tab[2][1] - tab[0][1] * tab[1][0] * tab[2][3] - tab[0][0] * tab[1][3] * tab[2][1] - tab[0][3] * tab[1][1] * tab[2][0];
		if ((w == 0) && (wx == 0) && (wy == 0) && (wz == 0))
		{
			cout << "Nieskonczenie wiele rozwiazan" << endl;
		}
		if ((w == 0) && (wx != 0 || wy != 0 || wz != 0))
		{
			cout << "Brak rozwiazan" << endl;
		}
		if (w != 0)
		{
			cout << "a2=" << (wx / w) << endl;
			cout << "a1=" << (wy / w) << endl;
			cout << "a0=" << (wz / w) << endl;
			cout << "Interpolowany wielomian :" << endl;
			cout << "W(x)=" << (wx / w) << "x^2" <<"+" << (wy / w) << "x" << (wz / w) << endl;
		}

	}
	plik.close();
}
void wielomianiowa()
{
	punkt *tab = new punkt[3];
	plik2.open("plik2.txt"), ios::out | ios::in;
	{
		int txtx, txty;
		for (int i = 0; i < 3; i++)
		{
			plik2 >> txtx >> txty;
			tab[i].x = txtx;
			tab[i].y = txty;
		}
		plik.open("plik.txt"), ios::out | ios::in;
		{
			for (int i = 0; i < 3; i++)
			{

				plik << (tab[i].x)*(tab[i].x) << " " << tab[i].x << " 1 " << tab[i].y << endl;
			}
		}
	}
	plik.close();
	plik2.close();
	wielomianowa2();
}

void Cramer()
{
	double tab[3][4];
	plik.open("plik.txt"), ios::out | ios::in;
	{
		double txtx, w = 0, wx, wy, wz;
		for (int j = 0; j < 3; j++)
		{
			for (int i = 0; i < 4; i++)
			{
				plik >> txtx;
				tab[j][i] = txtx;
			}
		}
		w = tab[0][0] * tab[1][1] * tab[2][2] + tab[0][1] * tab[1][2] * tab[2][0] + tab[0][2] * tab[1][0] * tab[2][1] - tab[0][1] * tab[1][0] * tab[2][2] - tab[0][0] * tab[1][2] * tab[2][1] - tab[0][2] * tab[1][1] * tab[2][0];
		wx = tab[0][3] * tab[1][1] * tab[2][2] + tab[0][1] * tab[1][2] * tab[2][3] + tab[0][2] * tab[1][3] * tab[2][1] - tab[0][1] * tab[1][3] * tab[2][2] - tab[0][3] * tab[1][2] * tab[2][1] - tab[0][2] * tab[1][1] * tab[2][3];
		wy = tab[0][0] * tab[1][3] * tab[2][2] + tab[0][3] * tab[1][2] * tab[2][0] + tab[0][2] * tab[1][0] * tab[2][3] - tab[0][3] * tab[1][0] * tab[2][2] - tab[0][0] * tab[1][2] * tab[2][3] - tab[0][2] * tab[1][3] * tab[2][0];
		wz = tab[0][0] * tab[1][1] * tab[2][3] + tab[0][1] * tab[1][3] * tab[2][0] + tab[0][3] * tab[1][0] * tab[2][1] - tab[0][1] * tab[1][0] * tab[2][3] - tab[0][0] * tab[1][3] * tab[2][1] - tab[0][3] * tab[1][1] * tab[2][0];

		for (int j = 0; j < 3; j++)
		{
			cout.setf(ios::showpos);
			cout << tab[j][0] << "x" << tab[j][1] << "y" << tab[j][2] << "z=" << tab[j][3] << endl;
		}
		cout << endl << endl;

		if ((w == 0) && (wx == 0) && (wy == 0) && (wz == 0))
		{
			cout << "Nieskonczenie wiele rozwiazan" << endl;
		}
		if ((w == 0) && (wx != 0 || wy != 0 || wz != 0))
		{
			cout << "Brak rozwiazan" << endl;
		}

		if (w != 0)
		{
			cout << "x=" << (wx / w) << endl;
			cout << "y=" << (wy / w) << endl;
			cout << "z=" << (wz / w) << endl;
		}

	}
	plik.close();
}

double function(double x,double max)
{
	double wsp,wartosc=0;
	plik2.open("plik2.txt"), ios::out | ios::in;
	{
		for (int i = 0; i < max+1; i++)
		{
			plik2 >> wsp;
			for(int j=0;j<i;j++)
			{
				wsp = wsp*x;
			}
			wartosc = wartosc + wsp;
		}
	}
	plik2.close();
	return wartosc;
}
double function2(double x, int max)
{
	double wartosc = 0, wsp;

	plik3.open("plik3.txt"), ios::out | ios::in;
	{
		for (int i = 0; i < max + 1; i++)
		{
			plik3 >> wsp;
			for (int j = i; j > 0; j--)
			{
				wsp = wsp*x;
			}
			wartosc = wartosc + wsp;
			wsp = 0;
		}
	}
	plik3.close();
	return wartosc;
}

double calki_prostokat()
{
	double max,xp,xk,n,dx,calka=0;
	cout << "Podaj stopien wielomianu" << endl;
	cin >> max;
	plik.open("plik.txt"), ios::out | ios::in;
	{
		plik >> xp;
		plik >> xk;
		plik >> n;
		dx = (xk - xp) / n;
	}
	cout << "Poczatek przedzialu: " << xp << ". Koniec przedzialu:  " << xk << ". Ilosc podzialow:  " << n << endl;
	for(int i=0;i<n;i++)
	{
		calka = calka + dx*function(xp+dx,max);
		xp = xp + dx;
	}
	plik.close();
	cout << "Calka wynosi: " << calka << endl;
	return calka;
}

double calki_trapez()
{
	double max, xp, xk, n, dx, calka=0;
	cout << "Podaj stopien wielomianu" << endl;
	cin >> max;
	plik.open("plik.txt"), ios::out | ios::in;
	{
		plik >> xp;
		plik >> xk;
		plik >> n;
		dx = (xk - xp) / n;
	}
	cout << "Poczatek przedzialu: " << xp << ". Koniec przedzialu:  " << xk << ". Ilosc podzialow:  " << n << endl;
	for (int i = 0;i<n;i++)
	{
		calka = calka + (function(xp , max)+function(xp+dx, max))*dx/2;
		xp = xp + dx;
	}
	plik.close();
	cout << "Calka wynosi: " << calka << endl;
	return calka;
}


double calki_simpson()
{
	double max, xp, xk, il_obsz, h, calka = 0;
	cout << "Podaj stopien wielomianu" << endl;
	cin >> max;
	plik.open("plik.txt"), ios::out | ios::in;
	{
		plik >> xp;
		plik >> xk;
		plik >> il_obsz;
		h = (xk - xp) / il_obsz;
	}
	cout << "Poczatek przedzialu: " << xp << ". Koniec przedzialu:  " << xk << ". Ilosc podzialow:  " << il_obsz << endl;
	for (int i = 0;i<(il_obsz/2);i++)
	{
		calka = calka + (function(xp, max) + 4*(function(xp + h, max))+ function(xp + h+h, max))*(h/3);
		xp = xp + h+h;
	}
	plik.close();
	cout << "Calka wynosi: " << calka << endl;
	return calka;
}


double calki_montecarlo()
{
	double max, xp, xk, il_losowan, calka = 0,random;
	cout << "Podaj stopien wielomianu" << endl;
	cin >> max;
	plik.open("plik.txt"), ios::out | ios::in;
	{
		plik >> xp;
		plik >> xk;
		plik >> il_losowan;
	}
	cout << "Poczatek przedzialu: " << xp << ". Koniec przedzialu:  " << xk << ". Ilosc losowan:  " << il_losowan << endl;

	for (int i = 0;i<il_losowan;i++)
	{
		random = (double)rand()/RAND_MAX*(xk - xp) + xp;
		calka = calka + function(random,max);
	}
	calka = (calka / il_losowan)*(xk-xp);
	plik.close();
	cout << "Calka wynosi: " << calka << endl;
	return calka;
}


double kwadratura_Gaussa()
{
	//punkty w pliku musza byc podane pokolei//
	punkt tab[4];
	double waga[] = { 1.0,1.0 };
	double punkt_Gaussa[] = { -0.5773502692,0.5773502692 };
	double N[2][2][4];
	double poch_ksi[2][4], poch_ni[2][4], matrix_det[2][2];
	plik.open("plik.txt"), ios::out | ios::in;
	{
		int txtx, txty;
		for (int i = 0; i < 4; i++)
		{
			plik >> txtx >> txty;
			tab[i].x = txtx;
			tab[i].y = txty;
			cout << "x:"<<txtx << " " <<"y:"<< txty << endl;
		}
	}
	plik.close();

	for (int j = 0; j < 2; j++)
	{
		for (int i = 0; i < 2; i++)
		{
			N[i][j][0] = 0.25*(1.0 - punkt_Gaussa[i])*(1.0 - punkt_Gaussa[j]);
			N[i][j][1] = 0.25*(1.0 + punkt_Gaussa[i])*(1.0 - punkt_Gaussa[j]);
			N[i][j][2] = 0.25*(1.0 + punkt_Gaussa[i])*(1.0 + punkt_Gaussa[j]);
			N[i][j][3] = 0.25*(1.0 - punkt_Gaussa[i])*(1.0 + punkt_Gaussa[j]);

			poch_ksi[j][0] = -0.25*(1.0 - punkt_Gaussa[j]);
			poch_ksi[j][1] = 0.25*(1.0 - punkt_Gaussa[j]);
			poch_ksi[j][2] = 0.25*(1.0 + punkt_Gaussa[j]);
			poch_ksi[j][3] = -0.25*(1.0 + punkt_Gaussa[j]);

			poch_ni[i][0] = -0.25*(1.0 - punkt_Gaussa[i]);
			poch_ni[i][1] = -0.25*(1.0 + punkt_Gaussa[i]);
			poch_ni[i][2] = 0.25*(1.0 + punkt_Gaussa[i]);
			poch_ni[i][3] = 0.25*(1.0 - punkt_Gaussa[i]);

		}
	}

	double dxksi, dyksi, dxni, dyni;
	for (int j = 0; j < 2; j++)
	{
		for (int i = 0; i < 2; i++)
		{
			dxksi = poch_ksi[j][0] * tab[0].x + poch_ksi[j][1] * tab[1].x + poch_ksi[j][2] * tab[2].x + poch_ksi[j][3] * tab[3].x;
			dyksi = poch_ksi[j][0] * tab[0].y + poch_ksi[j][1] * tab[1].y + poch_ksi[j][2] * tab[2].y + poch_ksi[j][3] * tab[3].y;
			dxni = poch_ni[i][0] * tab[0].x + poch_ni[i][1] * tab[1].x + poch_ni[i][2] * tab[2].x + poch_ni[i][3] * tab[3].x;
			dyni = poch_ni[i][0] * tab[0].y + poch_ni[i][1] * tab[1].y + poch_ni[i][2] * tab[2].y + poch_ni[i][3] * tab[3].y;
			matrix_det[i][j] = (dxksi*dyni) - (dxni*dyksi);
		}

	}
	double pole = 0;
	for (int j = 0; j < 2; j++)
	{
		for (int i = 0; i < 2; i++)
		{
			pole = pole + abs(matrix_det[i][j] * waga[i] * waga[j]);
		}

	}
	cout << "Pole powierzchni wynosi: " << pole << endl;
	return pole;
}

double nieliniowe_bisekcja()
{
	double a, b, zadana_dokladnosc = 0.0001;
	int max;
	cout << "podaj stopien wielomianu" << endl;
	cin >> max;

	plik.open("plik.txt"), ios::out | ios::in;
	{
		plik >> a;
		plik >> b;
	}
	plik.close();
	cout << "Szukanie rozwiazanie w przedziale [" << a << ";" << b << "]" << endl;
	cout << "Dokladnosc rozwiazania to " << zadana_dokladnosc << endl;
	double x1 = (a + b) / 2;

	while (abs(a - b)>zadana_dokladnosc)
	{
		if (function(x1, max) == 0)
		{
			return x1;
		}

		if (function(x1, max)*function(a, max) < 0)
			b = x1;
		if (function(x1, max)*function(b, max) < 0)
			a = x1;

		x1 = (a + b) / 2;
	}
	cout << "x=" << (a + b) / 2 << endl;
	return (a + b) / 2;

}
double nieliniowe_newtona_raphsona()
{
	double a, b, zadana_dokladnosc = 0.0001, x1 = 10;
	int max;
	cout << "podaj stopien wielomianu" << endl;
	cin >> max;

	plik.open("plik.txt"), ios::out | ios::in;
	{
		plik >> a;
		plik >> b;
	}
	plik.close();

	cout << "Szukanie rozwiazanie w przedziale [" << a << ";" << b << "]" << endl;
	cout << "Dokladnosc rozwiazania to " << zadana_dokladnosc << endl;

	while (true)
	{
		double x2 = x1 - (function(x1, max) / function2(x1, max - 1));
		if (abs(function(x1, max)) <= zadana_dokladnosc)
		{
			cout << "x=" << x1 << endl;
			return x1;
		}
		x1 = x2;
	}

}
void rk1() {
	double a, b, h, y0;
	int max, n, max2;
	cout << "podaj najwyzsza potege x" << endl;
	cin >> max;
	cout << "podaj najwyzsza potege y" << endl;
	cin >> max2;

	plik.open("plik.txt"), ios::out | ios::in;
	{
		plik >> a;
		plik >> b;
		plik >> h;
		plik >> y0;
	}
	cout << "x=<" << a << ";" << b << "> ." << "Krok h=" << h << endl;
	cout << "Warunek poczatkowy y(0)=" << y0 << endl;
	plik.close();
	n = (b - a) / h;
	punkt *tab = new punkt[n + 1];
	tab[0].x = a;
	tab[0].y = y0;
	for (int i = 1; i < n+1; i++)
	{
		tab[i].y = tab[i - 1].y + h*(function(tab[i - 1].x, max) + function2(tab[i - 1].y, max2));
		tab[i].x = tab[i - 1].x + h;

	}
	results.open("results.txt", ios::out | ios::app);
	{
		for (int i = 0; i < n+1; i++)
		{
			results << tab[i].x;
			results << "    ";
			results << tab[i].y << endl;
		}
		results << endl;
	}
	results.close();
}
void Huen()
{
	double a, b, h, y0;
	int max, n, max2;
	cout << "podaj najwyzsza potege x" << endl;
	cin >> max;
	cout << "podaj najwyzsza potege y" << endl;
	cin >> max2;

	plik.open("plik.txt"), ios::out | ios::in;
	{
		plik >> a;
		plik >> b;
		plik >> h;
		plik >> y0;
	}
	plik.close();
	n = (b - a) / h;
	punkt *tab = new punkt[n + 1];
	tab[0].x = a;
	tab[0].y = y0;
	for (int i = 1; i <= n; i++)
	{
		tab[i].y = tab[i - 1].y + (h / 2) * (function(tab[i - 1].x, max) + function2(tab[i - 1].y, max2) + function(tab[i - 1].x + h, max) + function2(tab[i - 1].y + h*function(tab[i - 1].x, max) + h*function2(tab[i - 1].y, max2), max2));
		tab[i].x = tab[i - 1].x + h;

	}
	results.open("results.txt", ios::out | ios::app);
	{
		for (int i = 0; i <= n; i++)
		{
			results << tab[i].x;
			results << " ";
			results << tab[i].y << endl;

		}
		results << endl << endl;
	}
	results.close();
}


void rk4()
{
	double a, b, h, y0, k1, k2, k3, k4;
	int max, n, max2;
	cout << "podaj najwyzsza potege x" << endl;
	cin >> max;
	cout << "podaj najwyzsza potege y" << endl;
	cin >> max2;

	plik.open("plik.txt"), ios::out | ios::in;
	{
		plik >> a;
		plik >> b;
		plik >> h;
		plik >> y0;
	}
	plik.close();
	n = (b - a) / h;
	punkt *tab = new punkt[n + 1];
	tab[0].x = a;
	tab[0].y = y0;
	for (int i = 1; i <= n; i++)
	{
		k1 = h*(function(tab[i - 1].x, max) + function2(tab[i - 1].y, max2));
		k2 = h*(function(tab[i - 1].x + 0.5*h, max) + function2(tab[i - 1].y + 0.5*k1, max2));
		k3 = h*(function(tab[i - 1].x + 0.5*h, max) + function2(tab[i - 1].y + 0.5*k2, max2));
		k4 = h*(function(tab[i - 1].x + h, max) + function2(tab[i - 1].y + k3, max2));
		tab[i].y = tab[i - 1].y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		tab[i].x = tab[i - 1].x + h;

	}
	results.open("results.txt", ios::out | ios::app);
	{
		for (int i = 0; i <= n; i++)
		{
			results << tab[i].x;
			results << " ";
			results << tab[i].y << endl;
		}
	}
	results.close();
}


void add_matrix()
{
	int max;
	cout << "Podaj wielkosc macierzy?" << endl;
	cin >> max;

	double **tab = new double*[max];
	double **tab2 = new double*[max];
	double **result = new double*[max];
	for (int i = 0; i < max; i++)
	{
		tab[i] = new double[max];
		tab2[i] = new double[max];
		result[i] = new double[max];
	}


		plik.open("plik.txt"), ios::out | ios::in;
		{
			for (int i = 0;i < max;i++)
			{
				for (int j = 0;j < max;j++)
				{
					plik >> tab[i][j];
				}
			}

			for (int i = 0;i < max;i++)
			{
				for (int j = 0;j < max;j++)
				{
					plik >> tab2[i][j];
				}
			}
			for (int i = 0;i < max;i++)
			{
				for (int j = 0;j < max;j++)
				{
					result[i][j] = tab[i][j] + tab2[i][j];
				}
			}

		}
		plik.close();
		for (int i = 0;i < max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << tab[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "+" << endl << endl;
		for (int i = 0;i < max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << tab2[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "=" << endl << endl;

		for (int i = 0;i < max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << result[i][j]<<" ";
			}
			cout << endl;
		}
	

}

void sub_matrix()
{
	int max;
	cout << "Podaj wielkosc macierzy?" << endl;
	cin >> max;

	double **tab = new double*[max];
	double **tab2 = new double*[max];
	double **result = new double*[max];
	for (int i = 0; i < max; i++)
	{
		tab[i] = new double[max];
		tab2[i] = new double[max];
		result[i] = new double[max];
	}

		plik.open("plik.txt"), ios::out | ios::in;
		{
			for (int i = 0;i < max;i++)
			{
				for (int j = 0;j < max;j++)
				{
					plik >> tab[i][j];
				}
			}

			for (int i = 0;i < max;i++)
			{
				for (int j = 0;j < max;j++)
				{
					plik >> tab2[i][j];
				}
			}
			for (int i = 0;i < max;i++)
			{
				for (int j = 0;j < max;j++)
				{
					result[i][j] = tab[i][j] - tab2[i][j];
				}
			}

		}
		for (int i = 0;i < max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << tab[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "-" << endl << endl;
		for (int i = 0;i < max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << tab2[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "=" << endl << endl;
		for (int i = 0;i < max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << result[i][j] << " ";
			}
			cout << endl;
		}
}
	


void multiply_matrix()
{
	int max;
	double suma = 0;
	cout << "Podaj wielkosc macierzy?" << endl;
	cin >> max;


	double **tab = new double*[max];
	double **tab2 = new double*[max];
	double **result = new double*[max];
	for (int i = 0; i < max; i++)
	{
		tab[i] = new double[max];
		tab2[i] = new double[max];
		result[i] = new double[max];
	}


	cout << endl;
		plik.open("plik.txt"), ios::out | ios::in;
		{
			for (int i = 0;i < max;i++)
			{
				for (int j = 0;j < max;j++)
				{
					plik >> tab[i][j];
				}
			}

			for (int i = 0;i <max;i++)
			{
				for (int j = 0;j < max;j++)
				{
					plik >> tab2[i][j];
				}
			}
			for (int i = 0; i < max; i++)
				for (int j = 0; j < max; j++)
				{
					suma = 0;
					for (int k = 0; k < max; k++)
						suma += tab[i][k] * tab2[k][j];
					result[i][j] = suma;
				}

		
		plik.close();
		for (int i = 0;i < max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << tab[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "*" << endl << endl;
		for (int i = 0;i < max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << tab2[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "=" << endl << endl;
		for (int i = 0;i <max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << result[i][j] << " ";
			}
			cout << endl;
		}
	}

}


void trans_matrix()
{
	int max;
	double suma = 0;
	cout << "Podaj wielkosc macierzy?" << endl;
	cin >> max;


	double **tab = new double*[max];
	double **result = new double*[max];
	for (int i = 0; i < max; i++)
	{
		tab[i] = new double[max];
		result[i] = new double[max];
	}


	cout << endl;
	plik.open("plik.txt"), ios::out | ios::in;
	{
		for (int i = 0;i < max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				plik >> tab[i][j];
			}
		}

		for (int i = 0; i < max; i++)
			for (int j = 0; j < max; j++)
			{
				result[i][j] = tab[j][i];
			}


		plik.close();
		for (int i = 0;i < max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << tab[i][j] << " ";
			}
			cout << endl;
		}

		cout << endl << "trans=>" << endl << endl;
		for (int i = 0;i <max;i++)
		{
			for (int j = 0;j < max;j++)
			{
				cout << result[i][j] << " ";
			}
			cout << endl;
		}
	}

}

int main()
{
	sub_matrix();
	_getch();
}