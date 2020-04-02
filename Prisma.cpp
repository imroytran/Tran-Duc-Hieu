#include <iostream>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <conio.h>
using namespace std;
#define PI 3.1415926535898

class Point {   // класс точки
	public:
	double x_, y_;     //координаты
	Point();
	void setP(double x, double y);
} A, B, C, D, E, F, G;
Point::Point(){
	x_ = 0.0;
	y_ = 0.0;
}
void Point::setP(double x, double y){
	x_ = x;
	y_ = y;
}

class Line {     // класс прямой линии
	public:
	double a_, b_, c_;       // уравнение прямой линии имеет вид Аx + By + C = 0
	Point X_, Y_;       // точки на прямой линии
	void Equa(Point X, Point Y);
	Line();
};
void Line::Equa(Point X, Point Y){          // Функция дает уравнение прямой линии когда изветстно две точки
		X_ = X;
		Y_ = Y;
		a_ = -(Y.y_ - X.y_);
		b_ = Y.x_ - X.x_;
		c_ = -(a_*X.x_ + b_*X.y_);
	}
Line::Line(){     // начальный лазерный луч
	a_ = 1.0; b_ = 0.0; c_ = -5.0;
}

class Prisma {
	double L_, l_;     // Параметры призмы
	double alpha_, beta_;
	public:
	Prisma();
	void setPr(double L, double l, double alpha, double beta);
	void Pris();
};
Prisma::Prisma(){
	L_ = 150.0;
	l_ = 15.0;
}
void Prisma::setPr(double L, double l, double alpha, double beta){
		L_ = L;
		l_ = l;
		alpha_ = alpha;
		beta_ = beta;
	}
void Prisma::Pris(){   // Функция дает все координаты всех основных точках треугольника
		A.setP(0.0, L_*tan(beta_));
		B.setP(-L_, 0.0);
		C.setP(L_, 0.0);
		E.setP(l_, 0.0);
		D.setP(-l_, 0.0);
		F.setP(0.0, l_*tan(alpha_));

	}

void EquaIn(Point X, Point Y, Line& t, Line k, double& gamma);
double Length(Point A, Point B);
bool Middle (Line k, Point G);
bool Same(Point G, Point H);
bool Control(Line t, Line k, Point G);
void PofInters(Line t, Line k);
bool ktra = false;
double Angle(double gamma, double beta);

int main()
{
	Line AB, AC, BC, EC, EF, DF;
	double gamma = 0.0;
	int j = 0;
	ofstream outfile;
	outfile.open("tran.txt");
	double L, l;
	cout << "Вводите значение L: ";
	cin >> L; cout << endl;
	cout <<"Вводите значение l: ";
	cin >> l;
	cout << endl;
	double tetamax;
	cout << "Вводите значение ограничения угла выхода: ";
	cin >> tetamax;
	Prisma PR1;
	cout << endl;
    //	double beta = (PI/2)*5.29/9;
     //	double alpha =(PI/2)*8.08/9;
    	for (double beta = (PI/2)*5.3/9; beta <= (PI/2)*8/9; beta += PI/1800){ // шаг угла beta
    		for(double alpha =(PI/2)*0/9 ; alpha < (PI/2)*8/9; alpha += PI/1800){ // шаг угла alpha

			PR1.setPr(L,l,alpha,beta);    // Ввод параметров призмы
			PR1.Pris();    // Создание призмы по основным точкам
			AB.Equa(A, B);  // уранение АВ
			BC.Equa(B, C);  // уранение BC
			AC.Equa(A, C);   // уранение AC
			EF.Equa(E, F);   // уранение EF
			DF.Equa(D, F);   // уранение DF
			Line t;    // Создание лазерного луча
			int m = 0;    // Количество переотражений
			PofInters(t, EF); // Нахождение первой точки пересечения луча с оптическим входом
			Point H;  // создание точки Н длч запоминания передыдущей точки пересечения
			H = G;   // Присвоивание координат точки пересечение в Н
			Line k, q;    //
			double r = asin(sin(alpha)/1.5);   // угол прелоления внутренней призмы
			double alpha1 = r - alpha;
			t.a_ = t.a_*cos(alpha1) - t.b_*sin(alpha1);        // Уранение луча после преломления
			t.b_ = t.a_*sin(alpha1) + t.b_*cos(alpha1);        //
			t.c_ = -(t.a_*G.x_ + t.b_*G.y_);                   //
                        k = AC;        // Присвоивание АС в линию к
                        PofInters(t, k);   // Поиск точки пересечения луча и АС
                        q = t;
                        m++;    // Увеличение количества переотражение на 1
                        EquaIn(H, G, t, AC, gamma);   // Поиск угла переотражения и уранения луча после переотражения
			Point T = G;   // самая высокая точка отражения
			ktra = false;
			
			while(Control(q, k, G) == false || (gamma >= asin(1/1.5)) || Angle(gamma, beta) > (tetamax/180)*PI || (G.y_ == 0)) {    // Условия продольжения цыкла
				H = G;
			        PofInters(t, DF);
				Point G1 = G;
				PofInters(t, EF);
				if ((Middle(DF, G1) == true && Same(G1,H) == false) || (Middle(EF, G) == true && Same(G, H) == false )){    // Если точка пересечения находится между двумя точками
					ktra = true;
					break;
				}
				for(int i = 0; i <3; i++){   // Проверка пересечения кождой стороны призмы
					switch(i){
						case 0:{
							k = AB;
							break;
						}
						case 1:{
							k = BC;
							break;
						}
                                                case 2:{
							k = AC;
							break;
						}
					}
					PofInters(t,k);
					if (Middle(k,G) == true && Same(G, H) == false){
						q = t;
						m++;
						if (G.y_ > T.y_) T.y_ = G.y_;
						EquaIn (H, G, t, k, gamma);
                                                break;
					}
				}
                                if (m > 7) break; // Изменение переотражений
				if (ktra == true) break;
			}
	     	        if (ktra == true || m > 7) continue;
	     	        if (A.y_ - T.y_ < 50*tan(beta)) continue;   // Условие что самая точка отражения не выше края сверху
	    	        double mm = (A.y_ - G.y_)/A.y_;
	    	        if (mm <= 0.25 || mm >= 0.75) continue;     // Лазерный луч только выходит в пределе [0,25H, 0,75H]
			cout << endl;
			cout << "Goc: " << (Angle(gamma, beta)/PI)*180 << endl;
			cout << "Diem cao nhat: " << T.y_ << endl;
			++j;
			outfile << (beta/PI)*180 << "          ";
			outfile << (alpha/PI)*180 << "          ";
			outfile << m << "         "; // Количество переотражений
			outfile << (fabs(asin(1.5*sin(gamma)) - (PI/2 - beta))/PI)*180 << endl;
			cout <<"m: " << m << endl;
		//	outfile << T.y_ << "           "; // самая высокая точка отражения
		//	outfile << G.y_ <<"         ";
		//	outfile << A.y_ << "         ";
			cout << "H: " << A.y_ << endl;
			cout << " beta: " << (beta/PI)*180;
			cout << " alpha: " << (alpha/PI)*180 << endl;
		//	if (G.x_ < 0) outfile << "лев" << endl;
		//	else outfile << "прав" << endl;
			cout <<"-------------------------------" << endl;
			ktra = false;
     	   	}
   	   	}

//	}
	cout << j;
	getch();
	return 0;
}

void EquaIn(Point X, Point Y, Line& t, Line k, double& gamma){     // Функция, дает уранение отражающего лича и уголь между ним и нормалью
	Point H0, H1;   // H0- точка проекции точки Н на нормали прямой линии k; H1 -  точка симетрии Н через нормаль прямой линии k в точки G
	H0.x_ = (k.a_*k.a_*X.x_ + k.b_*k.b_*Y.x_ + k.a_*k.b_*(X.y_ - Y.y_))/(k.a_*k.a_ + k.b_*k.b_);
	H0.y_ = (k.a_*X.x_ + k.b_*X.y_)/k.b_ - H0.x_*k.a_/k.b_;
	H1.x_ = 2*H0.x_ - X.x_;
	H1.y_ = 2*H0.y_ - X.y_;
	gamma = atan(Length(H1, H0)/Length(H0, Y));
	t.Equa(H1, G);
}

double Length(Point A, Point B){    // функция возвращается длины АВ
	return sqrt(pow(A.x_-B.x_,2)+pow(A.y_-B.y_,2));
}

bool Middle (Line k, Point G){    // Функция проверки пересечения G находится в нутрии призмы ли нет
        if ( k.X_.x_ > k.Y_.x_){
            if ( G.x_ > k.Y_.x_ && G.x_ < k.X_.x_ ) return true;
            else return false;
        }
        else{
            if ( G.x_ < k.Y_.x_ && G.x_ > k.X_.x_ ) return true;
            else return false;
		}
}
bool Same(Point G, Point H) {   // функция проверки совпадения точек
        if (fabs( H.x_ - G.x_) < 0.00001 && fabs(H.y_ - G.y_) < 0.00001) return true;
        else return false;

}

bool Control(Line t, Line k, Point G){  // функция условия приходящего лича находятся в выходящем диопазоне ли нет
	if (k.a_ == 0 || t.b_ == 0) return false;  // для исключения проверки горизонтальной стороны
	else{
    	if ((-t.c_/t.b_) > (k.a_*G.y_ - k.b_*G.x_)/k.a_) return true; // Точка пересечения лазерного лича и оси Oy выше чочки пересечения нормали через точку G прямой линии k и оси Oy
        else return false;
	}
}

void PofInters(Line t, Line k){   // Функция поиска координатов точки пересечения
      	double d, dx, dy;
		d = t.a_*k.b_ - k.a_*t.b_;
		dx = -t.c_*k.b_ + k.c_*t.b_;
		dy = -t.a_*k.c_ + k.a_*t.c_;
        if (d == 0) ktra = true ;
        else{
        G.x_ = dx/d;
        G.y_ =  dy/d;
        }
}
double Angle(double gamma, double beta){
	return fabs(asin(1.5*sin(gamma)) - (PI/2 - beta));
}





