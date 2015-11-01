clc;
disp( sprintf( '+--------------------------------------------------+' ) );
disp( sprintf( '|                                                  |  ' ) );
disp( sprintf( '|       EP - MAC0427 Programação Não Linear        |  ' ) );
disp( sprintf( '|              Marcos Kazuya Yamazaki              |  ' ) );
disp( sprintf( '|                                                  |  ' ) );
disp( sprintf( '|       NUSP: 7577622                              |  ' ) );
disp( sprintf( '|       Data: 01/Julho/2015                        |  ' ) );
disp( sprintf( '|       Professor: Walter Mascarenhas              |  ' ) );
disp( sprintf( '|                                                  |  ' ) );
disp( sprintf( '+--------------------------------------------------+\n' ) );

disp( sprintf( '   PROBLEMA: Minimização de uma função não linear' ) );
disp( sprintf( 'irretrita pelo método de Newnton regularizado\n' ) );

disp( sprintf( '   O programa aceita funções de até três variáveis,' ) );
disp( sprintf( 'sendo que elas tem que ser obrigatoriamente escritas' ) );
disp( sprintf( 'em funções das variavies x, y e z.' ) );
disp( sprintf( '   O programa pode nunca terminará caso o problema' ) );
disp( sprintf( 'é ilimitado, com o custo ótimo -infinito' ) );

clear;
syms x y z t;

disp( sprintf( '\nDigite a funcao que deseja minimizar' ) );
prompt = 'em função de x, y e z: ';
f = input(prompt);
display(f);

%   Variaveis que decidem o tamanho
% do passo a ser dado em cada etapa
%   Pode ser modificado, desde que 
% satisfaçam as suas restriçoes
% mencionadas abaixo.
alpha  =  0.5; % | 0 < alpha < 1
theta  =  0.5; % | 0 < theta < 1
lambda =  5.0; % |   lambda > 0
beta   =  5.0; % |   beta   > 0

disp( 'Digite os valores do primeiro ponto' );
prompt = 'x: '; X(1) = input(prompt); 
prompt = 'y: '; X(2) = input(prompt); 
prompt = 'z: '; X(3) = input(prompt);

% Calcular o valor da função Objetiva nesse ponto
valorObjetivo = subs(f            , x, X(1));
valorObjetivo = subs(valorObjetivo, y, X(2));
valorObjetivo = subs(valorObjetivo, z, X(3));

valorObjetivoAntigo = double(valorObjetivo);
disp(sprintf('Ponto de partida [%.3f, %.3f, %.3f]`', double(X(1)), double(X(2)),double(X(3))));
disp(sprintf('O valor da função neste ponto é %.4f\n', valorObjetivoAntigo));

%   Calcula o gradiente da função
% com as variáveis simbólicas.
G = [
    diff(f,x);
    diff(f,y);
    diff(f,z);
    ];

%   Calcula o Hessiano da função
% com as variaveis simbolicas.
H = [
    diff(G(1),x), diff(G(1),y), diff(G(1),z);
    diff(G(2),x), diff(G(2),y), diff(G(2),z); 
    diff(G(3),x), diff(G(3),y), diff(G(3),z); 
    ];

% Gradiente de f calculado no ponto X(1:3) atual;
g = G;
iteracao = 1;

while(1)
    disp(sprintf('------------------- Iteração %.0f -------------------\n', iteracao));
    iteracao = iteracao + 1;
    
    % Calcula o vetor gradiente no ponto atual
    for i = 1:3
        g(i) = subs(G(i), x, X(1));
        g(i) = subs(g(i), y, X(2));
        g(i) = subs(g(i), z, X(3));
    end
    disp(sprintf('Gradiente da  função: [%.3f, %.3f, %.3f]`', double(g(1)), double(g(2)),double(g(3))));
    % Caso se verifique que o gradiente é igual a zero
    % o ponto achado é um ponto de minimo local (ou máximo)
    if(norm(g) >= 0 && norm(g) <= 0.001 )
        disp(sprintf('\nPonto de mímino local: [%.3f, %.3f, %.3f]`', double(X(1)), double(X(2)),double(X(3))));
        disp(sprintf('Com valor da função: %.10f\n', double(valorObjetivo)));
        break;
    end

    % Calcula o Hessiano da função no ponto atual
    for i = 1:3
        for j = 1:3
            h(i,j) = subs(H(i,j), x, X(1));
            h(i,j) = subs(h(i,j), y, X(2));
            h(i,j) = subs(h(i,j), z, X(3));
        end
    end

    % Faz a fatoração de Cholesky para resolver o sistema lineares
    % afim de achar a direção d, na qual devemos procurar o proximo
    % ponto, que será o candidato a ser mínimo.
    [R,p] = chol(h);
    
    % Termos para verificar a condição do ângulo.
    % Coloquei esses valores apenas para entrar no loop seguinte.
    termo1 = 1;
    termo2 = 0;

    while(termo1 > termo2)
        %    Se na fatoração de Cholesky, o valor p retornado
        % for diferente de zero, isso quer dizer que a matriz
        % hessiana, h, não é positiva definida.
        %    Para ajustar isso, somamos por uma matriz identidade
        % multiplicada por um lambda. E tentamos novamente, e assim
        % consequentemente até a matriz h se tornar positiva definida.
        while (p ~= 0)
            h = h + lambda*eye(3);
            [R,p] = chol(h);
        end
        % Afim de achar a direção, d,
        % que se da pela formula: h*d = -g, onde h = R'*R
        % assim tempo R'*R*d = -g, onde R é uma matriz diagonal
        % inferior, que para se achar a sua inversa é muito mais 
        % rápida.
        % Fazemos: R'*r = -g, onde r = R*d =>
        %          r = -g*/R'
        %          d = r*/R
        r = vpa((-g'*inv(R'))');
        d = vpa((r'*inv(R))');
        
        % Verifica a condição do ângulo:
        termo1 = vpa(g'*d);
        termo2 = vpa(-(theta*norm(g)*norm(d)));
        
        p = 1; % Para o caso em que termo1 > termo2, pois dai
               % teriamos que entrar dentro do loop para calcular
               % a fatoração de Cholesky novamente.
        lambda = 10;
    end
    
    % Corrigir a direção
    if(norm(d) < beta*norm(g))
        d = (beta*norm(g)/norm(d))*d;
    end
    d = vpa(d);

    % Acha t, com a condicao de Armijo
    % Função composta, f(X + t*d)
    Xt = X(1) + t*d(1);
    Yt = X(2) + t*d(2);
    Zt = X(3) + t*d(3);
    
    valorObjetivoT = subs(f             , x, Xt);
    valorObjetivoT = subs(valorObjetivoT, y, Yt);
    valorObjetivoT = subs(valorObjetivoT, z, Zt);
    valorObjetivo  = subs(valorObjetivoT, t, 0 );
    outroLado = alpha*t*g'*d;
    
    % Acha um comprimeito do passo que satisfaça a condição de Armijo
    condicaoArmijo = valorObjetivo + outroLado - valorObjetivoT;
    comprimento    = sort(solve(condicaoArmijo));
    maiorSol       = size(comprimento);
    
    % Pega o passo de comprimento mais longo de valor absoluto
    % para evitar pegar o zero como o tamanho do passo
    
    passo = comprimento(maiorSol(1));
    if(passo <= 0.000001)
        disp(sprintf('\nA função é ilimitada'));
        disp(sprintf('Com o valor da função objetiva -infinito\n'));
        break;
    end
    
    % calcular o novo ponto!
    X(1) = X(1) + passo*d(1);
    X(2) = X(2) + passo*d(2);
    X(3) = X(3) + passo*d(3);

    % Calcular o valor da função Objetiva nesse novo ponto
    valorObjetivo = subs(f            , x, X(1));
    valorObjetivo = subs(valorObjetivo, y, X(2));
    valorObjetivo = subs(valorObjetivo, z, X(3));
    valorObjetivoAntigo = double(valorObjetivo);

    disp(sprintf('Passo e direção: %.3f*[%.2f, %.2f, %.2f]`\n', double(passo), double(d(1)), double(d(2)),double(d(3))));
    disp(sprintf('Novo ponto calculado: [%.3f, %.3f, %.3f]`', double(X(1)), double(X(2)),double(X(3))));
    disp(sprintf('O valor da função neste ponto é: %.9f\n', valorObjetivoAntigo));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EXEMPLO DE SAIDA ......................................................

%  Digite a funcao que deseja minimizar
%  em função de x, y e z: x^2 - 10*x + y^2 + 70*y
% 
%  f =
% 
%  x^2 - 10*x + y^2 + 70*y
% 
%  Digite os valores do primeiro ponto
%  x: 1
%  y: 1
%  z: 0
%  Ponto de partida [1.000, 1.000, 0.000]`
%  O valor da função neste ponto é 62.0000
%  
%  ------------------- Iteração 0 -------------------
%  
%  Gradiente da função: [-8.000, 72.000, 0.000]`
%  Passo e direção: 0.100*[40.00, -360.00, 0.00]`
%  
%  Novo ponto calculado: [5.000, -35.000, 0.000]`
%  O valor da função neste ponto é: -1250.000000000
%  
%  ------------------- Iteração 1 -------------------
%  
%  Gradiente da função: [0.000, 0.000, 0.000]`
%  
%  Ponto de mímino local: [5.000, -35.000, 0.000]`
%  Com valor da função: -1250.0000000000


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EXEMPLO DE SAIDA ......................................................
%  
%  Digite a funcao que deseja minimizar
%  em função de x, y e z: x^2 + 3*x*y + 5*y^2 + 2*x + 3*y^2
%  
%  f =
%   
%  x^2 + 3*x*y + 2*x + 8*y^2
%   
%  Digite os valores do primeiro ponto
%  x: 1
%  y: 1
%  z: 0
%  Ponto de partida [1.000, 1.000, 0.000]`
%  O valor da função neste ponto é 14.0000
%  
%  ------------------- Iteração 1 -------------------
%  
%  Gradiente da função: [7.000, 19.000, 0.000]`
%  Passo e direção: 0.012*[-33.15, -95.66, 0.00]`
%  
%  Novo ponto calculado: [0.595, -0.170, 0.000]`
%  O valor da função neste ponto é: 1.470548286
%  
%  ------------------- Iteração 2 -------------------
%  
%  Gradiente da função: [2.681, -0.929, 0.000]`
%  Passo e direção: 0.126*[-13.48, 4.42, 0.00]`
%  
%  Novo ponto calculado: [-1.109, 0.389, 0.000]`
%  O valor da função neste ponto é: -1.072238203
%  
%  ------------------- Iteração 3 -------------------
%  
%  Gradiente da função: [0.949, 2.894, 0.000]`
%  Passo e direção: 0.012*[-4.99, -14.39, 0.00]`
%  
%  Novo ponto calculado: [-1.170, 0.213, 0.000]`
%  O valor da função neste ponto é: -1.355731865
%  
%  ------------------- Iteração 4 -------------------
%  
%  Gradiente da função: [0.299, -0.104, 0.000]`
%  Passo e direção: 0.126*[-1.50, 0.49, 0.00]`
%  
%  Novo ponto calculado: [-1.360, 0.275, 0.000]`
%  O valor da função neste ponto é: -1.387338394
%  
%  ------------------- Iteração 5 -------------------
%  
%  Gradiente da função: [0.106, 0.323, 0.000]`
%  Passo e direção: 0.012*[-0.56, -1.60, 0.00]`
%  
%  Novo ponto calculado: [-1.367, 0.256, 0.000]`
%  O valor da função neste ponto é: -1.390862186
%  
%  ------------------- Iteração 6 -------------------
%  
%  Gradiente da função: [0.033, -0.012, 0.000]`
%  Passo e direção: 0.126*[-0.17, 0.05, 0.00]`
%  
%  Novo ponto calculado: [-1.388, 0.262, 0.000]`
%  O valor da função neste ponto é: -1.391255051
%  
%  ------------------- Iteração 7 -------------------
%  
%  Gradiente da função: [0.012, 0.036, 0.000]`
%  Passo e direção: 0.012*[-0.06, -0.18, 0.00]`
%  
%  Novo ponto calculado: [-1.389, 0.260, 0.000]`
%  O valor da função neste ponto é: -1.391298852
%  
%  ------------------- Iteração 8 -------------------
%  
%  Gradiente da função: [0.004, -0.001, 0.000]`
%  Passo e direção: 0.126*[-0.02, 0.01, 0.00]`
%  
%  Novo ponto calculado: [-1.391, 0.261, 0.000]`
%  O valor da função neste ponto é: -1.391303735
%  
%  ------------------- Iteração 9 -------------------
%  
%  Gradiente da função: [0.001, 0.004, 0.000]`
%  Passo e direção: 0.012*[-0.01, -0.02, 0.00]`
%  
%  Novo ponto calculado: [-1.391, 0.261, 0.000]`
%  O valor da função neste ponto é: -1.391304280
%  
%  ------------------- Iteração 10 -------------------
%  
%  Gradiente da função: [0.000, -0.000, 0.000]`
%  
%  Ponto de mímino local: [-1.391, 0.261, 0.000]`
%  Com valor da função: -1.3913042795


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EXEMPLO DE SAIDA ......................................................
%  
%  Digite a funcao que deseja minimizar
%  em função de x, y e z: x^3 + 3*x*y + 5*y^2 + 2*x + 3*y^2
%   
%  f =
%   
%  x^3 + 3*x*y + 2*x + 8*y^2
%   
%  Digite os valores do primeiro ponto
%  x: 1
%  y: 1
%  z: 0
%  Ponto de partida [1.000, 1.000, 0.000]`
%  O valor da função neste ponto é 14.0000
%  
%  ------------------- Iteração 1 -------------------
%  
%  Gradiente da função: [8.000, 19.000, 0.000]`
%  Passo e direção: 1.229*[-41.41, -94.39, 0.00]`
%  
%  Novo ponto calculado: [-49.901, -115.027, 0.000]`
%  O valor da função neste ponto é: -1291.862999607
%  
%  ------------------- Iteração 2 -------------------
%  
%  Gradiente da função: [7127.368, -1990.138, 0.000]`
%  
%  A função é ilimitada
%  Com o valor da função objetiva -infinito
