#ifndef TURBO_BASETYPES
#define TURBO_BASETYPES

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


inline double max(double a, double b) {
	return (a>b)?a:b;
};

inline double min(double a, double b) {
	return (a<b)?a:b;
};

inline int max(int a, int b) {
	return (a>b)?a:b;
};

inline int min(int a, int b) {
	return (a<b)?a:b;
};

inline std::pair<int, int> swap(std::pair<int,int>& pair)
{
	int buf = pair.first;
	pair.first = pair.second;
	pair.second = buf;
	return pair;
};

inline std::vector<double> abs(const std::vector<double>& a) {
	std::vector<double> res(5,0);
	for (int k = 0; k<a.size(); k++) res[k] = abs(a[k]);
	return res;
};

const double PI = 3.1459;

//Basic type definitions
//General exception
class Exception {
	std::string msg;
public:
	Exception(std::string m) {
		msg = m;
	};

	std::string ShowMessage() {
		std::cerr << msg.c_str() << "\n";
		return msg;
	};
};


//Vector type
class Vector {
public:
	double x;
	double y;
	double z;
	Vector(): x(0), y(0), z(0)		
	{
	};
	Vector(double _x, double _y, double _z) :
		x(_x),
		y(_y),
		z(_z)
		{
		};
	double mod() const
	{
		return sqrt(x*x+y*y+z*z);
	};

	double& operator[](int i) {
		if (i == 0) return x;
		if (i == 1) return y;
		if (i == 2) return z;
	};

	inline const Vector& operator+=(const Vector& a)
	{
		x += a.x;
		y += a.y;
		z += a.z;
		return *this;
	};
	inline const Vector& operator-=(const Vector& a)
	{
		x -= a.x;
		y -= a.y;
		z -= a.z;
		return *this;
	};
	inline const Vector& operator*=(const double& a)
	{
		x *= a;
		y *= a;
		z *= a;
		return *this;
	};	
	inline const Vector& operator/=(const double& a)
	{
		x /= a;
		y /= a;
		z /= a;
		return *this;
	};
};

inline double operator*(const Vector& a, const Vector& b)
{
	return (a.x*b.x + a.y*b.y + a.z*b.z);
};

inline Vector operator+(const Vector& a, const Vector& b)
{
	return Vector(a.x+b.x, a.y+b.y, a.z+b.z);
};

inline Vector operator-(const Vector& a, const Vector& b)
{
	return Vector(a.x-b.x, a.y-b.y, a.z-b.z);
};

inline Vector operator*(const double& a, const Vector& b)
{
	return Vector(a*b.x, a*b.y, a*b.z);
};

inline Vector operator*(const Vector& b, const double& a)
{
	return Vector(a*b.x, a*b.y, a*b.z);
};

inline Vector operator/(const Vector& b, const double& a)
{
	return Vector(b.x/a, b.y/a, b.z/a);
};

inline Vector operator-(const Vector& b)
{
	return Vector(-b.x, -b.y, -b.z);
};

inline Vector operator&(const Vector& a, const Vector &b)
{
	return Vector(a.y*b.z-b.y*a.z, b.x*a.z-a.x*b.z, a.x*b.y-a.y*b.x);
};

typedef std::vector<double> Row;
typedef std::vector<Row> Table;

//Matrix type
class Matrix {
private:
	int size_n;
	int size_m;
	Table table;
public:
	Matrix() {
	};

	Matrix(int n, int m) {
		size_n = n;
		size_m = m;
		table.resize(n , Row(m, 0));
	};

	Row& operator[](int i) {		
		return table[i];
	};	

	const Matrix& operator+=(Matrix& a) {
		for (int i = 0; i<size_n; i++) {
			for (int j = 0; j<size_m; j++) {
				table[i][j] += a[i][j];
			};
		};
		return *this;
	};

	const Matrix& operator*=(const double a) {
		for (int i = 0; i<size_n; i++) {
			for (int j = 0; j<size_m; j++) {
				table[i][j] *= a;
			};
		};
		return *this;
	};

	const Matrix& operator/=(const double a) {
		for (int i = 0; i<size_n; i++) {
			for (int j = 0; j<size_m; j++) {
				table[i][j] /= a;
			};
		};
		return *this;
	};
};

inline Matrix operator*(const double& a, const Matrix& b)
{
	Matrix res = b;
	res *= a;
	return res;
};

class RotationMatrix : public Matrix{
public:
	RotationMatrix() : Matrix(3,3) {		
		(*this)[0][0] = 1;
		(*this)[1][1] = 1;
		(*this)[2][2] = 1;
	};

	Vector operator*(Vector& x) {
		Vector res;
		for (int i = 0; i<3; i++) {
			for (int j = 0; j<3; j++) {
				res[i] += (*this)[i][j] * x[j];
			};
		};
		return res;
	};
};

//Class that represents shape
class Shape {
	
public:
	int marker;
	virtual bool isInside(Vector r) = 0;
};

class Box : public Shape {
private:
	Vector r_l;
	Vector r_r;
public:	
	Box(Vector rl, Vector rr) {
		if (rl.x > rr.x) throw Exception("Error in Box shape initialization");
		if (rl.y > rr.y) throw Exception("Error in Box shape initialization");
		if (rl.z > rr.z) throw Exception("Error in Box shape initialization");
		r_l = rl;
		r_r = rr;
	};

	bool isInside(Vector r) {
		if (!( (r_l.x <= r.x) && (r_r.x >= r.x))) return false;
		if (!( (r_l.y <= r.y) && (r_r.y >= r.y))) return false;
		if (!( (r_l.z <= r.z) && (r_r.z >= r.z))) return false;
		return true;
	};
};

class Circle : public Shape {
private:
	Vector r_c;
	double R;
public:
	Circle(Vector rc, double r) {
		r_c = rc;
		R = r;
	};

	bool isInside(Vector r) {
		return (r - r_c).mod() <= R;
	};

	void Translate(Vector dr) {
		r_c += dr;
	};
};

enum ShapeAdditionMode {
	Add,
	Subtract
};

class ComplexShape : public Shape {
private:
	std::vector<Shape*> shapes;
	std::vector<ShapeAdditionMode> modes;	
public:
	Vector r;
	Vector t;

	ComplexShape() {
		t = Vector(0,0,0);
		r = Vector(0,0,0);
	};

	void AddShape(Shape* s, ShapeAdditionMode m) {
		shapes.push_back(s);
		modes.push_back(m);
	};

	void Translate(Vector dr) {
		t += dr;
	};

	bool isInside(Vector _r) {
		bool isInside = false;
		Vector r = _r - t;
		for (int i = 0; i<shapes.size(); i++) {
			if (modes[i] == Add) {
				isInside |= shapes[i]->isInside(r);
			};
			if (modes[i] == Subtract) {
				isInside &= shapes[i]->isInside(r);
			};
		};
		return isInside;
	};
};

//Operators for vectors
template<class T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b) {
	if (a.size() != b.size()) throw Exception("Operation on different size vectors");
	std::vector<T> res(a.size());
	for (int i = 0; i<a.size(); i++) res[i] = a[i] + b[i];
	return res;
};

template<class T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b) {
	if (a.size() != b.size()) throw Exception("Operation on different size vectors");
	std::vector<T> res(a.size());
	for (int i = 0; i<a.size(); i++) res[i] = a[i] - b[i];
	return res;
};

template<class T>
std::vector<T> operator/=(std::vector<T>& a, const double c) {		
	for (int i = 0; i<a.size(); i++) a[i] /= c;
	return a;
};

template<class T>
std::vector<T> operator*(const std::vector<T>& a, const double c) {	
	std::vector<T> res(a.size());
	for (int i = 0; i<a.size(); i++) res[i] = a[i] * c;
	return res;
};

template<class T>
std::vector<T> operator/(const std::vector<T>& a, const double c) {	
	std::vector<T> res(a.size());
	for (int i = 0; i<a.size(); i++) res[i] = a[i] / c;
	return res;
};

//std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
//	if (a.size() != b.size()) throw Exception("Operation on different size vectors");
//	std::vector<double> res(a.size(), 0);
//	for (int i = 0; i<a.size(); i++) res[i] = a[i] + b[i];
//	return res;
//};

std::vector<double> operator*(const std::vector<double>& a, double c) {
	std::vector<double> res(a.size());
	for (int i = 0; i<a.size(); i++) res[i] = a[i] * c;
	return res;
};

std::vector<double> operator*(double c, const std::vector<double>& a) {		
	std::vector<double> res(a.size());
	for (int i = 0; i<a.size(); i++) res[i] = a[i] * c;
	return res;
};

std::vector<double> operator*=(std::vector<double>& a,double c) {		
	for (int i = 0; i<a.size(); i++) a[i] *= c;
	return a;
};

std::vector<double> operator+=(std::vector<double>& a, const std::vector<double>& b) {	
	if (a.size() != b.size()) throw Exception("Operation on different size vectors");
	for (int i = 0; i<a.size(); i++) a[i] += b[i];
	return a;
};

std::vector<double> operator-=(std::vector<double>& a, const std::vector<double>& b) {	
	if (a.size() != b.size()) throw Exception("Operation on different size vectors");
	for (int i = 0; i<a.size(); i++) a[i] -= b[i];
	return a;
};


//numerical function class
class function
{
public:
	int size;
	std::vector<double> x, value;
	function(int _size)
	{
		size = _size;
		x.resize(size);
		value.resize(size);
		for(int i=0; i<size; i++)
			value[i] = 0;
	};
	function(){};
	void ReadValue(char* filename); //just read values from file
	void ReadData1(char* filename);	//first x than values data arrangement in read file
	void ReadData2(char* filename);	//first values than x data arrangement in read file
	void WriteData(char* filename);	//write our function in a file
	void CreateUniformGrid(double l, double r, int N_nodes);
	double UniformSimpson();				//Integration by Simpson formula
	double TrapeziumIntegral();				//trapezium formula for non uniform grid
	function DerivateUniformGrid();			//derivation on uniform grid
	function Derivate();			//for nonuniform grid
	function DerivateSec();			//for nonuniform grid
	function BindToNewGridFO(std::vector<double>& grid); //interpolation of function f to the new grid
	function CutRight(double x_r);	//cut right part of function

	inline const function& operator+=(const function& a)
	{
		if(size!=a.size) throw Exception("functions have differrent size\n");
		for(int i=0; i<size; i++)
		{
			value[i] += a.value[i];
		};
		return *this;
	};
	inline const function& operator-=(const function& a)
	{
		if(size!=a.size) throw Exception("functions have differrent size\n");
		for(int i=0; i<size; i++)
		{
			value[i] -= a.value[i];
		};
		return *this;
	};
	inline const function& operator*=(const double& a)
	{
		for(int i=0; i<size; i++)
		{
			value[i] *= a;
		};
		return *this;
	};	
	inline const function& operator/=(const double& a)
	{
		for(int i=0; i<size; i++)
		{
			value[i] /= a;
		};
		return *this;
	};
	inline const function& operator+=(const double& a)
	{
		for(int i=0; i<size; i++)
		{
			value[i] += a;
		};
		return *this;
	};
	inline const function& operator-=(const double& a)
	{
		for(int i=0; i<size; i++)
		{
			value[i] -= a;
		};
		return *this;
	};
};

void function::ReadValue(char* filename)
{
	value.resize(0);
	x.resize(0);
	std::ifstream ifs(filename);
	double read;
	while(!ifs.eof())
	{
		ifs >> read;
		value.push_back(read);		//coordinate first
	};
	ifs.close();
	size = value.size();
	return;
};

void function::ReadData1(char* filename)
{
	value.resize(0);
	x.resize(0);
	std::ifstream ifs(filename);
	double read;
	while(!ifs.eof())
	{
		ifs >> read;
		x.push_back(read);		//coordinate first
		ifs >> read;
		value.push_back(read);
	};
	ifs.close();
	size = x.size();
	return;
};

void function::ReadData2(char* filename)
{
	value.resize(0);
	x.resize(0);
	std::ifstream ifs(filename);
	double read;
	while(!ifs.eof())
	{
		ifs >> read;
		value.push_back(read);		//value first
		ifs >> read;
		x.push_back(read);
	};
	size = x.size();
	ifs.close();
	return;
};

void function::CreateUniformGrid(double l, double r, int N_nodes)
{
	size = N_nodes;
	x.resize(size);
	for(int i=0; i<size; i++)
	{
		x[i] = l + i*(r - l)/(N_nodes - 1);
	};
};

void function::WriteData(char* filename)
{
	std::ofstream ofs(filename);
	int i;
	for(i=0; i<x.size()-1; i++)
	{
		ofs << x[i] << " " << value[i] << "\n";
	};
	ofs << x[i] << " " << value[i];
	ofs.close();
	return;
};

double function::UniformSimpson()
{
	int N = x.size()-1;
	if(N%2!=0) return TrapeziumIntegral();
	double A = 0;
	double B = 0;
	for(int i=1; i<=N-1; i += 2) A += 4*value[i];
	for(int i=2; i<=N-2; i += 2) B += 2*value[i];
	double res = value[0] + 4*A + 2*B + value[N];
	res /= 3.0*N;
	return res;
};

double function::TrapeziumIntegral()
{
	double I = 0;
	for (int i=1; i<x.size(); i++)
		I += 0.5*(value[i] + value[i-1])*(x[i] - x[i-1]);
	return I;
};

function function::DerivateUniformGrid()
{
	function Diff(value.size());
	Diff.value[0] = (value[1] - value[0])/(x[1] - x[0]);
	Diff.x[0] = x[0];
	for(int i=1; i<Diff.size-1; i++)
	{
		Diff.value[i] = (value[i+1] - value[i-1])/(x[i+1] - x[i-1]);
		Diff.x[i] = x[i];
	};
	Diff.value[Diff.size-1] = (value[Diff.size-1] - value[Diff.size-2])/(x[Diff.size-1] - x[Diff.size-2]);
	Diff.x[Diff.size-1] = x[Diff.size-1];
	return Diff;
};

function function::Derivate()
{
	function Diff(value.size());
	Diff.value[0] = (value[1] - value[0])/(x[1] - x[0]);			//left point
	Diff.x[0] = x[0];
	//second order derivatives approximation
	for(int i=1; i<Diff.size-1; i++)
	{
		double r = x[i+1] - x[i];
		double l = x[i] - x[i-1];
		Diff.value[i] = (l*l*value[i+1] - r*r*value[i-1] + (r*r - l*l)*value[i])/(l*r*(l + r));
		Diff.x[i] = x[i];
	};
	//right point
	Diff.value[Diff.size-1] = (value[Diff.size-1] - value[Diff.size-2])/(x[Diff.size-1] - x[Diff.size-2]);
	Diff.x[Diff.size-1] = x[Diff.size-1];
	return Diff;
};

function function::DerivateSec()
{
	function Diff(value.size());
	if(Diff.size<3) throw Exception("There is not enough function size to compute second derivative!\n");
	//left point
	double h = x[2] - x[0];
	double h1 = x[1] - x[0];
	Diff.value[0] = 2.0*(value[1]*h - value[2]*h1 - value[0]*(h - h1))/(h*h1*(h1 - h));
	Diff.x[0] = x[0];
	//second order derivatives approximation
	for(int i=1; i<Diff.size-1; i++)
	{
		double r = x[i+1] - x[i];
		double l = x[i] - x[i-1];
		Diff.value[i] = 2.0*(value[i+1]*l + value[i-1]*r - value[i]*(l+r))/(l*r*(l+r));
		Diff.x[i] = x[i];
	};
	//right point
	h = x[Diff.size-1] - x[Diff.size-3];
	h1 = x[Diff.size-1] - x[Diff.size-2];
	Diff.value[Diff.size-1] = 2.0*(value[Diff.size-2]*h - value[Diff.size-3]*h1 - value[Diff.size-1]*(h-h1));
	Diff.value[Diff.size-1] /= h*h1*(h1-h);
	Diff.x[Diff.size-1] = x[Diff.size-1];
	return Diff;
};

function function::BindToNewGridFO(std::vector<double>& grid)
{
	if((this->size==0)||(grid.size()==0)) throw new Exception("Function size is zero\n");

	//new function
	function res(0);
	double during_point = grid[0];
	int left = 0;		//index of left interval border
	double r = (this->x[1]) - (this->x[0]);
	double r_l, r_r;
	if(during_point < (this->x[left])) throw new Exception("Can't interpolate function to the new grid\n");
	for(int i=0; i<grid.size(); i++){
		during_point = grid[i];
		if((this->x[left + 1]) >= during_point){			
			r_l = (during_point - this->x[left])/r;
			r_r = (this->x[left + 1] - during_point)/r;
			res.value.push_back(r_l*this->value[left + 1] + r_r*this->value[left]);
			res.x.push_back(during_point);		
			continue;
		}
		left++;
		r = (this->x[left + 1]) - (this->x[left]);
		i--;
	};
	res.size = res.x.size();
	return res;
};

function function::CutRight(double x_r)
{
	function res;
	for(int i=0; (this->x[i])<=x_r; i++)
	{
		res.value.push_back(this->value[i]);
		res.x.push_back(this->x[i]);
	};
	return res;
};

function operator*(function a, double c) {		
	for (int i = 0; i<a.size; i++) a.value[i] *= c;
	return a;
};

function operator*(double c, function a) {		
	for (int i = 0; i<a.size; i++) a.value[i] *= c;
	return a;
};

function operator+(function a, double c) {		
	for (int i = 0; i<a.size; i++) a.value[i] += c;
	return a;
};

function operator+(double c, function a) {		
	for (int i = 0; i<a.size; i++) a.value[i] += c;
	return a;
};

function operator-(function a, double c) {		
	for (int i = 0; i<a.size; i++) a.value[i] -= c;
	return a;
};

function operator-(double c, function a) {		
	for (int i = 0; i<a.size; i++) a.value[i] = c - a.value[i];
	return a;
};

function operator/(function a, double c) {		
	for (int i = 0; i<a.size; i++) a.value[i] /= c;
	return a;
};

function operator+(function a, function& b) {		
	if(a.size!=b.size) throw Exception("functions have different sizes\n");
	for (int i = 0; i<a.size; i++) a.value[i] += b.value[i];
	return a;
};

function operator-(function a, function& b) {		
	if(a.size!=b.size) throw Exception("functions have different sizes\n");
	for (int i = 0; i<a.size; i++) a.value[i] -= b.value[i];
	return a;
};

//double type with global index
struct RealValue
{
public:
	int GlobalIndex;
	double value;

	RealValue(int _GlobalIndex) : GlobalIndex(_GlobalIndex) {};
	RealValue(){};
};

//tensor type with global index
struct TensorValue
{
public:
	int GlobalIndex;
	Matrix value;

	TensorValue(int _GlobalIndex) : GlobalIndex(_GlobalIndex) {};
	TensorValue(){};
};






#endif