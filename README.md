# Burgers-equation-FDM-solution
Straight method againts flow:
$$
\frac{u_n^{i+1} - u_n^i}{\Delta t} + u \frac{u_n^i - u_{n-1}^i}{\Delta x} = \nu \frac{u_{i+1}^n - 2u_i^n + u_{i-1}^n}{\Delta x^2}
$$

![Явная схема против потока](https://github.com/user-attachments/assets/25fb6e19-e5bc-45fd-87b4-defc7a17904d)

Lax method:
$$
\frac{u_n^{i+1}-0.5\ (u_{n+1}^i-u_{n-1}^i)}{dt}+u\frac{u_{n+1}^i-u_{n-1}^i}{2dx}=\nu(\frac{u_{i+1}-2u_i-u_{i-1}}{dx^2})
$$
![Схема Лакса](https://github.com/user-attachments/assets/320db663-7cb0-41f6-9cad-d955bd5288eb)

MacCormack’s method:
$$
\frac{{\hat{u}}_n^{i+1}-u_n^i}{dt}+u\frac{u_{n+1}^i-u_n^i}{dx}=\nu(\frac{u_{i+1}-2u_i-u_{i-1}}{dx^2})
\frac{{\bar{\hat{u}}}_j^{i+1}-u_n^i}{dt}+u\frac{{\hat{u}}_n^{i+1}-{\hat{u}}_{n-1}^{i+1}}{dx}=\nu(\frac{{\hat{u}}_{n+1}^i-2*{\hat{u}}_n^i-{\hat{u}}_{n-1}^i}{dx^2})
u_n^{i+1}=0.5\ ({\hat{u}}_n^i+{\bar{\hat{u}}}_j^{i+1})
$$
![Схема Маккормака](https://github.com/user-attachments/assets/3f76ddb1-c997-41ed-88c7-5a0bb263577f)

Two steps Lax-Wendroff method:
$$
\frac{u_{n+0.5}^\ast-0.5\left(u_{n+1}^i-u_{n-1}^i\right)}{dt/2}+u\frac{u_{n+1}^i-u_{n-1}^i}{dx}=\nu(\frac{u_{i+1}-2u_i-u_{i-1}}{dx^2})
\frac{u_n^{i+1}-u_n^i}{dt}+u\frac{u_{n+1/2}^\ast-u_{n-1/2}^\ast}{dx}=0
$$
![Схема Лакса-Вендроффа](https://github.com/user-attachments/assets/b1feadb9-68bc-4a12-a466-b6efd2dd2f23)

Conclusions: solution of parabolic type PDE by using FDM methods requires to respect CFL, correct choose viscosity coefficient, because wrong choose could raise dissipation and oscillation of the solution. All graph solution can find in folder.
