/**********************************************
* Este programa es la heuristica para TSP *
***********************************************/

#define MAXCOUNT 10 + 1
#include <iostream>
#include <string>
#include <cassert>
//using namespace std;      

// int main(){
// 	// // Imprimir Mensajes
// 	// std::cout "Hello World\n";

// 	// // Operaciones Basicas
// 	// float x;
// 	// x = 5.0/2.0;
// 	// std::cout << x;

// 	// // Vectores
// 	// float zip[5];
// 	// zip[0] = 1.0/3.0;
// 	// std::cout << zip[0] << "\n";

// 	// // Strings (Usar paquete <string>)
// 	// int length_name;
// 	// char first_char;
// 	// std::string my_name;
// 	// my_name = "Alfredo";
// 	// first_char = my_name.at(0);
// 	// length_name = my_name.length();
// 	// std::cout << first_char << "\n";
// 	// std::cout << length_name << "\n";

// 	// // Lectura de inputs de usuario
// 	int price, number_on_hand;
// 	std::cout << "Escribe los siguientes valores.\nPrice: ";
// 	std::cin >> price;
// 	std::cout << "Number on hand: ";
// 	std::cin >> number_on_hand;
// 	return (0);
// }

// int  height;   /* the height of the triangle */
// int  width;    /* the width of the triangle */
// int  area;     /* area of the triangle (computed) */

// int main(  )
// {
//     std::cout << "Enter width height? ";
//     std::cin >> width >> height;
//     area = (width * height) / 2;
//     std::cout << "The area is " << area << '\n';
//     return (0);
// }

// enum {
//       LOOKUP =1, /* default - looking rather than defining. */
//       VERB,
//       ADJ,
//       ADV,
//       NOUN,
//       PREP,
//       PRON,
//       CONJ
// };

// void my_function(int *a_ptr, int *b_ptr){
// 	int temp;
// 	temp = *a_ptr;
// 	*a_ptr = *b_ptr;
// 	*b_ptr = temp;
// }

// int main(){
//     int a, b;
//     a = 2;
//     b = 3;
//     my_function(&a, &b);
//     printf("%i \n", b);
//     printf("%i \n", LOOKUP);
//     return (0);
// }