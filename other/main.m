cor = true;
w = 576;
h = 768;
background = zeros([w h 3]);
background = background_func(cor);
cores= ['y','m', 'c', 'r', 'g', 'b', 'w', 'k'];
corPessoa = 'g';
corCarro = 'r';
corUndefined = 'b';

while true

close all, clc
if (cor)
   disp('O vнdeo vai ser visualizado a cores.');
else
    disp('O vнdeo vai ser visualizado a preto e branco.');
end
disp('1: Trocar entre cores e preto e branco');
disp('2: Ver vнdeo original');
disp('3: Ver background');
disp('4: Ver vнdeo com detecзгo de objectos');
disp('5: Ver vнdeo com detecзгo de objectos com altura e largura');
disp('6: Ver vнdeo com detecзгo de objectos com бrea');
disp('7: Alterar cor da detecзгo de pessoas');
disp('8: Alterar cor da detecзгo de carros');
disp('9: Alterar cor da detecзгo de desconhecidos');
disp('10: Sair');
opcao = input('Escolha a opзгo: ');
switch(opcao)
    case 1
        if(cor)
            disp('Trocou de cor para preto e branco. O background vai ser calculado novamente.');
                    background = zeros([w h]);
        else
            disp('Trocou de preto e branco para cor. O background vai ser calculado novamente. ');
                    background = zeros([w h 3]);
        end
        cor = not(cor);
        disp('Carregue espaзo para calcular');
        pause;
        background = background_func(cor);
    case 2
        showvideo(cor);
    case 3
        imshow(background), colormap gray;
    case 4
        find_objects(cor, background, corPessoa, corCarro, corUndefined, false, false);
    case 5
        find_objects(cor, background, corPessoa, corCarro, corUndefined, true, false);
    case 6
        find_objects(cor, background, corPessoa, corCarro, corUndefined, false, true);
    case 7
        disp('Qual a cor que quer para a pessoa?');
        disp('1) Amarelo');
        disp('2) Magenta');
        disp('3) Cian');
        disp('4) Vermelho');
        disp('5) Verde');
        disp('6) Azul');
        disp('7) Branco');
        disp('8) Preto');
        corescolhida = input('Escolha a cor: ');
        corPessoa = cores(corescolhida);
    case 8
        disp('Qual a cor que quer para o carro?');
        disp('1) Amarelo');
        disp('2) Magenta');
        disp('3) Cian');
        disp('4) Vermelho');
        disp('5) Verde');
        disp('6) Azul');
        disp('7) Branco');
        disp('8) Preto');
         corescolhida = input('Escolha a cor: ');
        corCarro = cores(corescolhida);
    case 9
        disp('Qual a cor que quer para o desconhecido?');
        disp('1) Amarelo');
        disp('2) Magenta');
        disp('3) Cian');
        disp('4) Vermelho');
        disp('5) Verde');
        disp('6) Azul');
        disp('7) Branco');
        disp('8) Preto');
        corescolhida = input('Escolha a cor: ');
        corUndefined =  cores(corescolhida);
    case 10
        return;
    otherwise
        disp('Opcгo invбlida');
    end
    fprintf('Carregue espaзo para continuar\n');
    pause;
end
