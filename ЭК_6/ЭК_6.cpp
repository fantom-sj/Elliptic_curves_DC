#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <type_traits>
#include <assert.h>
using namespace std;

template<class T1, class T2>
inline T1 mod(T1 a, T2 b) {
	assert(b > 0);
	T1 ret = a % b;
	if constexpr (std::is_unsigned_v<T1>)
	{
		return ret;
	}
	else {
		return (ret >= 0) ? (ret) : (ret + b);
	}
}

struct Point {
	long x;
	long y;

	Point() = default;
	Point(long point_x, long point_y) {
		x = point_x;
		y = point_y;
	}

	friend std::ostream& operator<< (std::ostream& out, const Point& point) {
		out << "(" << point.x << ", " << point.y << ")";
		return out;
	}
	friend bool operator== (const Point& A, const Point& B) {
		if ((A.x == B.x) && (A.y == B.y))
			return true;
		else
			return false;
	}
};

class ElipticCurv {
private:
	long a, b;
	long p;
	Point P = Point(1, 1);

public:
	ElipticCurv() = default;
	ElipticCurv(long p_a, long p_b) {
		a = p_a;
		b = p_b;
		p = 0;
	}
	ElipticCurv(long p_a, long p_b, long mod_p) {
		if (Prime(mod_p) == 1) {
			a = p_a;
			b = p_b;
			p = mod_p;
		}
		else {
			cout << "Инициализировать элептическую кривую с данным значеним модуля невозможно" << endl;
		}
	}
	Point sum_point(Point A, Point B) {
		if (A == B)
			return sum_point(A);
		else {
			Point N;
			long chisl = mod(B.y - A.y, p);
			long znam = mod(B.x - A.x, p);
			long inv_znam = mod(ext_evclid(znam, p), p);
			long l = mod(chisl * inv_znam, p);

			N.x = mod(powm(l, 2, p) - A.x - B.x, p);
			N.y = mod((l * (A.x - N.x)) - A.y, p);
			return N;
		}
	}

	Point sum_point(Point A) {
		Point N;
		long chisl = mod(3 * powm(A.x, 2, p) + a, p);
		long znam = mod(2 * A.y, p);
		long inv_znam = mod(ext_evclid(znam, p), p);
		long l = mod(chisl * inv_znam, p);

		N.x = mod(powm(l, 2, p) - (2 * A.x), p);
		N.y = mod((l * (A.x - N.x)) - A.y, p);
		return N;
	}
	Point point_mult_n(Point A, long n) {			//Умножение числа на точку в элептических кривых
		vector<long> bits = expand_two(n);

		vector<Point> mass_point;
		for (int i = 0; i < bits.size(); i++) {
			mass_point.push_back(pow_point(A, bits[i]));
		}

		Point sum = mass_point[0];
		for (int i = 1; i < mass_point.size(); i++) {
			sum = sum_point(sum, mass_point[i]);
		}
		return sum;
	}

	long ext_evclid(long x, long y) {
		long a, b, q, r, u1, u2, t;

		a = x;
		b = y;
		u1 = 1;
		u2 = 0;

		while (b != 0) {
			q = a / b;
			r = a % b;
			a = b; b = r;
			t = u2;
			u2 = u1 - q * u2;
			u1 = t;
		}
		return u1;
	}


private:
	long Prime(long a) {			//Проверка простоты числа
		long i1, i2, i3, i4, i5, i6, i7, i8, bound;

		if (a == 0 || a == 1) return 0;
		if (a == 2 || a == 3 || a == 5 || a == 7 || a == 11 || a == 13 || a == 17 || a == 19 || a == 23 || a == 29) return 1;
		if (a % 2 == 0 || a % 3 == 0 || a % 5 == 0 || a % 7 == 0 || a % 11 == 0 || a % 13 == 0 || a % 17 == 0 || a % 19 == 0 || a % 23 == 0 || a % 29 == 0) return 0;

		bound = sqrt(a);
		i1 = 31; i2 = 37; i3 = 41; i4 = 43; i5 = 47; i6 = 49; i7 = 53; i8 = 59;

		while (i8 <= bound && a % i1 && a % i2 && a % i3 && a % i4 && a % i5 && a % i6 && a % i7 && a % i8) {
			i1 += 30; i2 += 30; i3 += 30; i4 += 30; i5 += 30; i6 += 30; i7 += 30; i8 += 30;
		}

		if (i8 <= bound ||
			i1 <= bound && a % i1 == 0 ||
			i2 <= bound && a % i2 == 0 ||
			i3 <= bound && a % i3 == 0 ||
			i4 <= bound && a % i4 == 0 ||
			i5 <= bound && a % i5 == 0 ||
			i6 <= bound && a % i6 == 0 ||
			i7 <= bound && a % i7 == 0)
			return 0;

		return 1;
	}
	vector<long> expand_two(long n) {			//Разложение n по степеням числа 2
		vector<long> deg;
		while (n) {
			long bin_power = n & ~(n - 1), pos_bit = 0;
			while (bin_power) {
				pos_bit++;
				bin_power = bin_power >> 1;
			}
			deg.push_back(pos_bit - 1);
			n &= n - 1;
		}
		return deg;
	}
	Point pow_point(Point A, long stepen) {			//Возведение точки в степень
		Point N;
		if (stepen == 0)
			return A;
		else if (stepen == 1)
			return sum_point(A);
		else {
			N = sum_point(pow_point(A, stepen - 1));
		}
		return N;
	}
	long powm(long a, long n, long m) {
		if (n == 0)
			return 1 % m;
		if (n % 2 == 1)
			return (powm(a, n - 1, m) * a) % m;
		else
			return powm((a * a) % m, n / 2, m);
	}
};

Point cod_ascii_point[159] = {
	Point(33,355), Point(33,396), Point(33,74), Point(34,677), Point(36,87), Point(36,664), Point(39,171),
	Point(39,580), Point(43,224), Point(43,527), Point(44,366), Point(44,385), Point(45,31), Point(45,720),
	Point(47,349), Point(47,402), Point(48,49), Point(48,702), Point(49,183), Point(49,568), Point(53,277),
	Point(53,474), Point(56,332), Point(56,419), Point(58,139), Point(58,612), Point(59,365), Point(59,386),
	Point(61,129), Point(61,622), Point(62,372), Point(62,379), Point(66,379), Point(66,552), Point(67,84),
	Point(67,667), Point(69,241), Point(69,510), Point(70,195), Point(70,556), Point(72,254), Point(72,497),
	Point(73,72), Point(73,679), Point(74,170), Point(74,581), Point(75,318), Point(75,433), Point(78,271),
	Point(78,480), Point(79,111), Point(79,640), Point(80,318), Point(80,433), Point(82,270), Point(82,481),
	Point(83,373), Point(83,378), Point(85,35), Point(85,716), Point(86,25), Point(86,726), Point(90,21),
	Point(90,730), Point(93,267), Point(93,484), Point(98,338), Point(98,413), Point(99,295), Point(99,456),
	Point(100,364), Point(10,387), Point(102,267), Point(102,484), Point(105,369), Point(105,382), Point(106,24),
	Point(106,727), Point(108,247), Point(108,504), Point(109,200), Point(109,551), Point(110,129), Point(110,622),
	Point(114,144), Point(114,607), Point(115,242), Point(115,509), Point(116,92), Point(116,659), Point(120,147),
	Point(120,604), Point(125,292), Point(125,459), Point(126,33), Point(189,297), Point(189,454), Point(192,32),
	Point(192,719), Point(194,205), Point(194,546), Point(197,145), Point(197,606), Point(198,224), Point(198,527),
	Point(200,30), Point(200,721), Point(203,324), Point(203,427), Point(205,372), Point(205,379), Point(206,106),
	Point(206,645), Point(209,82), Point(209,669), Point(210,31), Point(210,720), Point(215,247), Point(215,504),
	Point(218,150), Point(218,601), Point(221,138), Point(221,613), Point(226,9), Point(226,742), Point(227,299),
	Point(227,452), Point(228,271), Point(228,480), Point(229,151), Point(229,600), Point(234,164), Point(234,587),
	Point(235,19), Point(235,732), Point(236,39), Point(236,712), Point(237,297), Point(237,454), Point(238,175),
	Point(238,576), Point(240,309), Point(240,442), Point(243,87), Point(243,664), Point(247,266), Point(247,485),
	Point(249,183), Point(249,568), Point(250,14), Point(250,737), Point(251,245), Point(251,506), Point(253,211),
	Point(253,540), Point(256,121), Point(256,630), Point(257,293), Point(257,458)
};
char ascii_table[159] = {
	' ', '!', '"', '#', '$', '%', '&', '\'', '(', ')', ' * ', ' + ', ',', ' - ', '.', ' / ',
	'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', '; ', ' < ', ' = ', ' > ', ' ? ',
	'@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q',
	'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a', 'b', 'c',
	'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u',
	'v', 'w', 'x', 'y', 'z', '{', '|', '}', '~', 'А', 'Б', 'В', 'Г', 'Д', 'Е', 'Ж', 'З', 'И',
	'Й', 'К', 'Л', 'М', 'Н', 'О', 'П', 'Р', 'С', 'Т', 'У', 'Ф', 'Х', 'Ц', 'Ч', 'Ш', 'Щ', 'Ъ',
	'Ы', 'Ь', 'Э', 'Ю', 'Я', 'а', 'б', 'в', 'г', 'д', 'е', 'ж', 'з', 'и', 'й', 'к', 'л', 'м',
	'н', 'о', 'п', 'р', 'с', 'т', 'у', 'ф', 'х', 'ц', 'ч', 'ш', 'щ', 'ъ', 'ы', 'ь', 'э', 'ю', 'я'
};

class ECP_gen {
private:
	ElipticCurv E = ElipticCurv(-1, 1, 751);

public:
	int r = 0,
		s = 0;
	Point Q;

	ECP_gen() = default;
	ECP_gen(int r_, int s_, Point Q_) {
		r = r_;
		s = s_;
		Q = Q_;
	}
	ECP_gen(int e, int d, int k, Point G, int n) {
		Point kG;

		do {
			kG = E.point_mult_n(G, k);
			// cout << kG << endl;
			r = mod(kG.x, n);
			// cout << r << endl;
		} while (r == 0);

		int Z = E.ext_evclid(k, n);
		// cout << Z << endl;
		
		s = Z * (e + d * r);
		s = mod(s, n);

		Q = E.point_mult_n(G, d);
	};
};

bool verefi_ecp(int e, ECP_gen ecp, int n, Point G) {
	ElipticCurv E = ElipticCurv(-1, 1, 751);
	if (!(1 <= ecp.r <= n - 1) ||
		!(1 <= ecp.s <= n - 1))
		return false;
	
	// cout << ecp.r << endl;
	int V = mod(E.ext_evclid(ecp.s, n), n);
	// cout << V << endl;
	int u1 = mod(e * V, n);
	// cout << u1 << endl;
	int u2 = mod(ecp.r * V, n);
	// cout << u2 << endl;
	Point x1 = E.point_mult_n(G, u1);
	// cout << "x1 " << x1 << endl;
	Point x2 = E.point_mult_n(ecp.Q, u2);
	// cout << "x2 " << x2 << endl;
	Point X = E.sum_point(x1, x2);
	// cout << X << endl;

	// cout << mod(X.x, n) << endl;
	if (mod(X.x, n) == ecp.r)
		return true;
	else
		return false;
}

int main() {
	setlocale(LC_ALL, "Russian");

	/*-------------------Создание подписи---------------------*/
	Point G = Point(416, 55);
	int n = 13;
	int e = 3,
		d = 9,
		k = 6;

	ECP_gen ecp = ECP_gen(e, d, k, G, n);
	cout << "Была сгенерирована следующая цифровая подпись: (" << ecp.r << ", " << ecp.s << ")"
		<< " С открытым ключём: " << ecp.Q << endl;
	

	
	/*-------------------Проверка подписи---------------------*/
	// G = Point(562, 89);
	// n = 13;
	// e = 5;
	// ecp.r = 3;
	// ecp.s = 7;
	// ecp.Q = Point(455, 368);
	bool ver_ecp = verefi_ecp(e, ecp, n, G);
	cout << "Была проверена следующая цифровая подпись: (" << ecp.r << ", " << ecp.s << ")"
		<< " С открытым ключём: " << ecp.Q << endl;
	if (ver_ecp)
		cout << "Подпись является верной!" << endl;
	else
		cout << "Подпись была откланена!" << endl;


	cout << endl;
	system("pause");
	return 0;
}

/*int e[10] = { 6, 6, 3, 8, 9, 6, 10, 6, 11, 10 };
int d[10] = { 7, 5, 2, 12, 6, 10, 9, 10, 5, 5 };
int k[10] = { 11, 6, 8, 8, 6, 2, 11, 7, 7, 11 };
ECP_gen ecp[10];

for (int i = 0; i < 10; i++) {
	ecp[i] = ECP_gen(e[i], d[i], k[i], G, n);
	cout << "Была сгенерирована следующая цифровая подпись: (" << ecp[i].r << ", " << ecp[i].s << ")"
		<< "С открытым ключём: " << ecp[i].Q << endl;
	bool ver_ecp = verefi_ecp(e[i], ecp[i], n, G);
	if (ver_ecp)
		cout << "Подпись является верной!" << endl;
	else
		cout << "Подпись была откланена!" << endl;
	cout << endl;
}*/