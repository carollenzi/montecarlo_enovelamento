#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

double energia(int N_e, double J_e[210], int *rede_e, int *cadeia_e);

int main()
{

    FILE *frede, *fenergia, *ftamanho, *fvar; // energia e tamanho ponta a ponta da cadeia em funcao do tempo (até 5e5 mct), guardar a rede a cada 1e5
    // media da energia e do tamanho em funcao da temperatura e variacao da energia em funcao da temperatura : 20e5 mct em cada T
    frede = fopen("rede_temperatura.txt", "w");
    fenergia = fopen("media_energia_temperatura.txt", "w");
    ftamanho = fopen("media_tamanho_temperatura.txt", "w");
    fvar = fopen("variacao_energia_temperatura.txt", "w");

    int N = 15;
    double T = 10;

    int ca[15] = {13, 5, 16, 10, 1, 13, 2, 4, 20, 18, 5, 9, 3, 6, 14};

    // cadeia de N aminoacidos (proteina)
    int *cadeia = (int *)malloc(N * sizeof(int));
    srand(time(NULL));
    for (int i = 0; i < N; i++)
    {
        // cada posicao tem um aminoacido aleatorio dentre os 20 possiveis
        // *(cadeia + i) = rand() % 20 + 1;
        *(cadeia + i) = ca[i];
    }

    // matriz da rede / arranjo espacial quadrado para colocar a proteina
    // cada posicao da rede vai guardar um indice da cadeia
    int *rede = (int *)malloc((N * N) * sizeof(int));

    // inicializar a rede (-1 default, pra diferenciar dos valores dos indices da cadeia de aminoacidos (de 0 a N-1))
    memset(rede, -1, (N * N) * sizeof(int));

    // comecar a proteina esticada no meio (vertical) da rede
    int inicio_cadeia = N * (int)floor(N / 2);
    int cont_c = 0;

    for (int i = inicio_cadeia; i < inicio_cadeia + N; i++)
    {
        *(rede + i) = cont_c;
        cont_c++;
    }

    // calcular energia J relativa aos pares de aminoacidos

    static double J[210] = {-2.805016, -2.930263, -3.100754, -3.479467, -2.391199, -3.376290, -2.440084, -3.118355, -3.830022, -2.074231, -2.685897, -3.495036, -3.583785, -2.509137, -2.621861, -2.126574, -3.649594, -3.650929, -2.231168, -3.663477, -3.507747, -3.415482, -3.290485, -3.760605, -2.576387, -2.553644, -2.841172, -2.271862, -3.457017, -3.030782, -2.487337, -2.262033, -3.961046, -3.588092, -3.741500, -2.352245, -2.964381, -2.181584, -3.470600, -2.794403, -2.255816, -2.156497, -2.289439, -3.839601, -2.665634, -2.911300, -3.966175, -2.315228, -2.562230, -2.197343, -3.978705, -2.069977, -3.612825, -3.269190, -3.830582, -2.189212, -3.822834, -2.671754, -2.461075, -3.279851, -3.702537, -2.948412, -3.541885, -3.663582, -2.536503, -3.283385, -2.015827, -3.500885, -3.464969, -3.486427, -2.295288, -3.720785, -3.642924, -2.584728, -3.560386, -2.308558, -3.496028, -3.526561, -2.623786, -2.058258, -3.723904, -2.602491, -2.128234, -3.336729, -3.871681, -3.958816, -3.525942, -3.694515, -2.630570, -3.987016, -2.974367, -2.333107, -2.935428, -2.516251, -3.996689, -3.471931, -3.799636, -2.012517, -2.972816, -3.264605, -3.498944, -3.268104, -2.985390, -3.141868, -3.852832, -2.545775, -3.450426, -3.348860, -2.072336, -2.074211, -3.407117, -3.796240, -2.676702, -3.535352, -3.132970, -2.548384, -3.494168, -2.658911, -2.242899, -2.124738, -2.645927, -3.217266, -2.457845, -3.581355, -3.733517, -2.454534, -3.053287, -3.533153, -2.467051, -2.026103, -2.797758, -3.965995, -3.294207, -3.783148, -3.107863, -3.147039, -2.328923, -2.558288, -2.495899, -2.401260, -2.632500, -3.903016, -2.197500, -3.309202, -3.438368, -3.330470, -3.857586, -2.932535, -3.989381, -2.100485, -3.057273, -2.635309, -3.317751, -3.515118, -2.216664, -3.051268, -3.969653, -3.269951, -2.584421, -2.436704, -3.296054, -3.382179, -2.402699, -2.590261, -3.165327, -3.510562, -3.737300, -3.494250, -2.068850, -2.233199, -3.895510, -2.701350, -2.136215, -2.093010, -2.010552, -3.574583, -3.423480, -3.868138, -2.507118, -3.412862, -3.968623, -3.564391, -2.048170, -3.286374, -3.079510, -2.264834, -2.337641, -3.049162, -3.534785, -2.922062, -3.485866, -2.830839, -2.304241, -3.888565, -3.421100, -3.469567, -3.399127, -3.158400, -2.963818, -3.467977, -3.391598, -2.859328, -2.169327, -3.527813, -2.952338, -2.179879, -3.102396, -2.375818, -2.048017, -3.609514};
    // double min = -4.0, max = -2.0;  // teste com numeros aleatorios entre -4 e -2 como no livro
    /*
    for (int i = 1; i <= 210; i++)
    {
        J[i] = min + (rand() / (double)RAND_MAX) * (max - min);
    }
    */

    // comecar enovelamento

    int mct = 0; // monte carlo time
    int t_max = (int)20e5;
    int aminoacido;  // indice da cadeia (um aminoacido para mover)
    int posicao;     // posicao do aminoacido na rede
    int vizinhos[8]; // vetor de vizinhos da posicao sorteada na rede (0: cima, 1: baixo, 2: dir, 3: esq, 4: cima dir, 5: cima esq, 6: baixo dir, 7: baixo esq)
    int v;           // indice no vetor de vizinhos
    bool l1, l2, mov_aceito, v_existe;
    double Ei, Ef, deltaE, prob_dE, rd, E_media = 0.0, E2_media = 0.0, varE, E_mct;
    int pos_inicio, pos_fim, tam; // posicoes das pontas da cadeia para calcular o tamanho
    double tam_media = 0.0;

    while (T > 0) // temperature sweep
    {
        mct = 0;

        while (mct <= t_max)
        {

            if (mct % (int)5e5 == 0)
            {
                fprintf(frede, "\nrede em %d mct, T = %f\n", mct, T);
                for (int i = 0; i < N * N; i++)
                {
                    fprintf(frede, "%d  ", *(rede + i));
                    if ((i + 1) % N == 0)
                        fprintf(frede, "\n");
                }
            }

            aminoacido = rand() % N;

            for (int i = 0; i < (N * N); i++)
            {
                if (*(rede + i) == aminoacido)
                {
                    posicao = i;
                    break;
                }
            }

            // se vizinhos[i] == -1 a posicao sorteada nao tem o vizinho i
            for (int i = 0; i < 8; i++)
            {
                vizinhos[i] = -1;
            }

            if (posicao - N >= 0) // cima
            {
                vizinhos[0] = posicao - N;
            }

            if (posicao + N < N * N) // baixo
            {
                vizinhos[1] = posicao + N;
            }

            if (posicao + 1 < N * N && (abs((posicao + 1) % N - (posicao % N))) == 1) // direita
            {
                vizinhos[2] = posicao + 1;
            }

            if (posicao - 1 >= 0 && (abs((posicao - 1) % N - (posicao % N))) == 1) // esquerda
            {
                vizinhos[3] = posicao - 1;
            }

            if (posicao - N >= 0 && (posicao - N + 1 % N) != 0) // cima e direta
            {
                vizinhos[4] = posicao - N + 1;
            }

            if (posicao - N - 1 >= 0 && (posicao - N - 1 % N) != (N - 1)) // cima e esquerda
            {
                vizinhos[5] = posicao - N - 1;
            }

            if (posicao + N + 1 < N * N && (posicao + N + 1) != 0) // baixo e direta
            {
                vizinhos[6] = posicao + N + 1;
            }

            if (posicao + N < N * N && (posicao + N - 1 % N) != (N - 1)) // baixo e esquerda
            {
                vizinhos[7] = posicao + N - 1;
            }

            // energia na configuracao inicial (antes do movimento)
            Ei = energia(N, J, rede, cadeia);

            do
            {
                // sortear uma posicao vizinha para ver se e possivel trocar o aminoacido de lugar
                v = rand() % 8;

                // l1 e l2: verficar se ligacoes com aminoacidos a esq e a dir nao serao quebradas (esticadas ou comprimidas) no movimento
                l1 = false;
                l2 = false;
                mov_aceito = false;

                // verificar se existe o vizinho sorteado
                if (vizinhos[v] != -1)
                {
                    v_existe = true;

                    // verificar se a posicao do vizinho na rede está vazia
                    if (*(rede + vizinhos[v]) == -1)
                    {
                        // verificar se é possivel mover o aminoacido para vizinhos[v] sem comprimir ou esticar as ligacoes entre aminoacidos
                        // os aminoacidos ligados ao que talvez se mova ainda estao na vizinhanca 4-conectada da nova posicao?

                        // analisar o aminoacido ligado a esquerda - l1
                        if (aminoacido - 1 >= 0)
                        {
                            if (*(rede + vizinhos[v] + N) == aminoacido - 1 || *(rede + vizinhos[v] - N) == aminoacido - 1 ||
                                *(rede + vizinhos[v] + 1) == aminoacido - 1 || *(rede + vizinhos[v] - 1) == aminoacido - 1)
                            {
                                l1 = true;
                            }
                        }
                        else
                            l1 = true; // aminoacido a esquerda nao existe (ponta da cadeia) => apenas l2 deve ser decisivo para o movimento

                        // analisar o aminoacido ligado a direita - l2
                        if (aminoacido + 1 < N)
                        {
                            if (*(rede + vizinhos[v] + N) == aminoacido + 1 || *(rede + vizinhos[v] - N) == aminoacido + 1 ||
                                *(rede + vizinhos[v] + 1) == aminoacido + 1 || *(rede + vizinhos[v] - 1) == aminoacido + 1)
                            {
                                l2 = true;
                            }
                        }
                        else
                            l2 = true; // aminoacido a direita nao existe (ponta da cadeia) => apenas l1 deve ser decisivo para o movimento

                        if (l1 && l2)
                        {
                            // mover aminoacido e "esvaziar" a posicao original
                            *(rede + vizinhos[v]) = aminoacido;
                            *(rede + posicao) = -1;

                            // calcular energia da proteina na nova configuracao
                            Ef = energia(N, J, rede, cadeia);

                            deltaE = Ef - Ei;

                            if (deltaE >= 0)
                            {
                                prob_dE = exp(-deltaE / T); // kB = 1
                                rd = (double)rand() / (double)RAND_MAX;

                                if (rd < prob_dE)
                                    mov_aceito = true;
                                else
                                    mov_aceito = false;
                            }
                            else
                                mov_aceito = true;

                            if (!mov_aceito) // desfazer o movimento
                            {
                                *(rede + vizinhos[v]) = -1;
                                *(rede + posicao) = aminoacido;
                            }
                        }
                    }
                }
                else
                    v_existe = false;

                E_mct = mov_aceito ? Ef : Ei;
                E_media += E_mct;
                E2_media += E_mct * E_mct;

                // energia em funcao do tempo
                // fprintf(fenergia, "%d   %f\n", mct, E_mct);

                // tamanho ponta a ponta em funcao do tempo (posicao inicio cadeia - posicao fim cadeia)

                for (int i = 0; i < (N * N); i++)
                {
                    if (*(rede + i) == 0)
                    {
                        pos_inicio = i;
                    }

                    if (*(rede + i) == (N - 1))
                    {
                        pos_fim = i;
                    }
                }

                tam = abs(pos_inicio - pos_fim);
                tam_media += tam;

                // fprintf(ftamanho, "%d   %d\n", mct, tam);

                mct++; // contar um tempo a cada sorteio (considerar tambem quando nao e possivel fazer o movimento)

            } while (!v_existe);
        }

        E_media /= mct;
        E2_media /= mct;
        varE = mct / (mct - 1) * (E2_media - (E_media * E_media));
        tam_media /= mct;

        fprintf(fenergia, "%f   %f\n", T, E_media);
        fprintf(ftamanho, "%f   %f\n", T, tam_media);
        fprintf(fvar, "%f   %f\n", T, varE);

        T -= 0.5;
    }

    fclose(frede);
    fclose(fenergia);
    fclose(ftamanho);
    fclose(fvar);

    return 0;
}

double energia(int N_e, double J_e[210], int *rede_e, int *cadeia_e)
{
    // E = sum(delta(i, j) * J(i, j))
    // delta(i, j) = 1 se aminoacidos i e j sao vizinhos na rede mas nao estao ligados (nao sao adjacentes na cadeia)
    // delta(i, j) = 0 caso contrario

    int pos_i, a, ij, ji;
    double E = 0.0;

    // percorrer pares de aminoacidos da cadeia
    for (int i = 0; i < N_e; i++)
    {
        // encontrar aminoacido i na rede
        for (int k = 0; k < N_e * N_e; k++)
        {
            if (*(rede_e + k) == i)
            {
                pos_i = k;
                break;
            }
        }

        for (int j = 0; j < N_e; j++)
        {
            if (j != i) // se j == i é o mesmo aminoacido
            {
                if (j != i + 1 && j != i - 1) // i e j nao estao ligados
                {
                    // ver se i e j sao vizinhos na rede
                    if (*(rede_e + pos_i + N_e) == j || *(rede_e + pos_i - N_e) == j || *(rede_e + pos_i + 1) == j || *(rede_e + pos_i - 1) == j)
                    {
                        // Aij = J(i + ((2n-j)(j-1)/2)), i>=j
                        ij = *(cadeia_e + i) + ((40 - *(cadeia_e + j)) * (*(cadeia_e + j) - 1) / 2);
                        ji = *(cadeia_e + j) + ((40 - *(cadeia_e + i)) * (*(cadeia_e + i) - 1) / 2);

                        a = *(cadeia_e + i) >= *(cadeia_e + j) ? ij : ji;
                        E += J_e[a]; // calcular a energia relativa aos aminoacidos nas posicoes i e j
                    }
                }
            }
        }
    }

    return E;
}