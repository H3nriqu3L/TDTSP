#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include <chrono>
#include <queue>
#include <random>


using namespace std;

struct PairWithValue {
    int key;
    int value;

    // Construtor
    PairWithValue(int k, int v) : key(k), value(v) {}

    bool operator<(const PairWithValue& other) const {
        return value > other.value;  // Ordenação em ordem crescente pelo valor associado
    }
};

struct NewRota
{
    vector<int> rota;
    int i,j;
    int custo;
    NewRota(const vector<int>& r, int x, int y, int c) : rota(r), i(x), j(y), custo(c) {}
    NewRota(): NewRota({}, 0, 0,0) {}
};

// Declaracoes funcoes
void buscaLocalBestImp(const vector<vector<int>> &dist, vector<int> &rota, const vector<vector<int>> &multas);
void grasp(const vector<vector<int>> &dist, vector<int> &rota, const vector<vector<int>> &multas, int maxSegundos, float alp);
vector<int> geraRotaAleatoria(int n);
void buscaLocalBestImpTempo(const vector<vector<int>> &dist, vector<int> &rota, const vector<vector<int>> &multas, int maxSegundos);

int calcDistEuclides(float xa, float ya, float xb,float yb){
    return round(sqrt( pow((xb-xa),2) + pow((yb-ya),2)));
}

int calcDistEntreCidades(pair<float,float> c1, pair<float,float> c2){
    return calcDistEuclides(c1.first, c1.second, c2.first, c2.second);
}

int calcResultado(const vector<int> &lista,const vector<vector<int>> &dist, const vector<vector<int>> &multas ){ // Lista é nossa rota, em que o elemento lista[0] esta ligado em lista[1]
    int sumResult = 0;
    
    // Itera a lista (rota)
    for(int i=0; i<lista.size(); i++){
        //Caso i seja a ultima cidade somamos o custo da ultima cidade ate a primeira para fechar o ciclo
        //cout << sumResult << " ";
        if(i==(lista.size()-1)){
            //sumResult+=calcDistEntreCidades(coord[lista[i]], coord[lista[0]]);
            sumResult+=dist[lista[i]][lista[0]];
            sumResult+= multas[lista[i]][i];
            break;
        }

        //sumResult+=calcDistEntreCidades(coord[lista[i]], coord[lista[i+1]]);
        sumResult+=dist[lista[i]][lista[i+1]];
        sumResult+= multas[lista[i]][i];
    }

    return sumResult;
}

int calcDifMulta(int i,const vector<int> &rota, const vector<vector<int>> &multas){     //Calcula a multa de uma rota, sendo i a iteracao que a primeira cidade foi visitada
    int sumMulta = 0;
    for(int k=0; k<rota.size(); k++){
        sumMulta+=multas[rota[k]][i+k];
        //cout << "A multa cobrada pela cidade " << rota[k] << " por ser a " << i+k << " cidade ée de " << multas[rota[k]][i+k] << endl;
    }

    return sumMulta;
}

void printRota(const vector<int> &rota){
    for(int i=0; i<rota.size(); i++){
        cout << rota[i]+1 << " ";
    }
    cout << endl;
}

vector<vector<int>> geraMatrizDist(const vector<pair<float,float>> &coord){ //gera matriz[0][1] = distancia cidade 0 ate 1 ou vice versa
    int tam = coord.size();
    vector<vector<int>> dist(tam, vector<int>(tam, 0.0));

    for (int i = 0; i < tam; i++) {
        for (int j = i + 1; j < tam; j++) {
            int distancia = calcDistEntreCidades(coord[i],coord[j]);
            dist[i][j] = distancia;
            dist[j][i] = distancia;
        }
    }

    return dist;
}

int encontrarProximaCidade(const vector<vector<int>>& distancias, vector<bool>& visitado, int cidadeAtual) {    //Encontra a cidade mais proxima dado uma cidade (nao considera multa)
    int proximaCidade = -1;
    int distanciaMinima = numeric_limits<int>::max();

    for (int i = 0; i < distancias[cidadeAtual].size(); i++) {
        if (!visitado[i] && distancias[cidadeAtual][i] < distanciaMinima) {
            proximaCidade = i;
            distanciaMinima = distancias[cidadeAtual][i];
        }
    }

    return proximaCidade;
}

priority_queue<PairWithValue> encontraRLCProximaCidade(const vector<vector<int>>& distancias, vector<bool>& visitado,const vector<vector<int>> &multas, int cidadeAtual, int iteracao) {    //Encontra a cidade mais proxima dado uma cidade (nao considera multa)
    vector<int> cidades;
    priority_queue<PairWithValue> pq;
    int n = visitado.size();
    vector<int> sequencia = geraRotaAleatoria(n); //Iremos comparar as cidades de ordem aleatoria pois em grafos grandes nao conseguimos comparar toda a vizinhanca
    int max_iteracoes = 600, cont=0;


    for (int i = 0; i < n; i++) {
        if(cont>max_iteracoes && !pq.empty()) break; //Em grafos muito grandes fica inviavel olhar toda a vizinhanca, entao olharemos somente uma parte
        if (!visitado[sequencia[i]] ) {
            pq.push(PairWithValue(sequencia[i], (distancias[cidadeAtual][sequencia[i]]+multas[sequencia[i]][iteracao])));
        }
        cont++;
    }

    return pq;
}


vector<int> geraRotaGuloso(const vector<vector<int>> &dist, int size){
    int n = size;
    vector<bool> visitado(n, false);
    vector<int> rota(n); //Um vetor que guardará a ordem que devemos visitar as cidades
    rota[0] = 0;   //Primeira será a 0
    visitado[0] = true;


    for(int i=1; i<n; i++){
        int cidadeAtual = rota[i - 1];
        int proximaCidade = encontrarProximaCidade(dist, visitado, cidadeAtual);
        visitado[proximaCidade] = true;
        rota[i] = proximaCidade;

    }

    return rota;

}

vector<int> geraRotaGuloso2(const vector<vector<int>> &dist, int size ){
    int n = size;
    vector<bool> visitado(n, false);
    vector<int> rota; //Um vetor que guardará a ordem que devemos visitar as cidades

    rota.push_back(0);
    visitado[0] = true;


    // Acrescenta a cidade mais próxima da cidade 0
    int proximaCidade = encontrarProximaCidade(dist, visitado, 0);
    rota.push_back(proximaCidade);
    visitado[proximaCidade] = true;

    for(int i=2; i<n; i++){
        int cidadeComeco = rota[0];
        int cidadeFim = rota[rota.size()-1];

        int proximaCidadeComeco = encontrarProximaCidade(dist, visitado, cidadeComeco);
        int proximaCidadeFim = encontrarProximaCidade(dist, visitado, cidadeFim);

        if(dist[cidadeComeco][proximaCidadeComeco]<dist[cidadeFim][proximaCidadeFim]){
            rota.insert(rota.begin(), proximaCidadeComeco);
            visitado[proximaCidadeComeco]=true;
        }else{
            rota.push_back(proximaCidadeFim);
            visitado[proximaCidadeFim]=true;
        }
    }

    return rota;
}

vector<int> geraRotaSequencial(int size){
    //Gera sequencia para testes da P0
    int n = size;
    vector<int> teste(n);
    for(int i=0; i<n; i++){
        teste[i] = i;
    }
    return teste;
}

vector<int> geraRotaAleatoria(int n) {
    //Gera sequencia para testes da P0
    vector<int> rota(n);
    rota[0] = 0;
    for(int i=1; i<n; i++){
        rota[i] = i;
    }
    random_device rd;
    mt19937 rng(rd());
    shuffle(rota.begin()+1, rota.end(), rng);

    return rota;
}

vector<int> geraSequenciaAleatoria(int inicio, int n){
    //Gera sequencia para testes da P0
    vector<int> rota(n);
    for(int j=0; j<inicio; j++){
        rota[j]=0;
    }
    for(int i=inicio; i<n; i++){
        rota[i] = i;
    }
    random_device rd;
    mt19937 rng(rd());
    shuffle(rota.begin(), rota.end(), rng);

    return rota;
}

vector<NewRota> gerarMovimentosVizinhos(const vector<int>& rota, const int custoRota, const vector<vector<int>> &dist, const vector<vector<int>> &multas,  const vector<vector<int>> &listaTabu) {
    vector<NewRota> vizinhos;
    int n = rota.size();
    //int custo = calcResultado(rota, dist, multas);
    int custo = custoRota;

    //Otimizacao para grafos muito grandes
    bool otimizar=false;
    int contador=0;
    float chanceJump=0;
    if(n>300){
        otimizar=true;
        chanceJump = 0.7;

    }
    auto startTime = std::chrono::steady_clock::now();
    

    for (int i = 1; i < n-1; i++) {
        for (int j = i + 1; j < n; j++) {
            if(listaTabu[i][j]!=0) continue;


            NewRota vizinho;
            vizinho.rota = rota;
            reverse(vizinho.rota.begin() + i, vizinho.rota.begin() +j + 1);
            //vizinhos.push_back(vizinho);
            vizinho.i = i; vizinho.j = j;


            vector<int> multaRotaPreInversao, multaRotaPosInversao;
            multaRotaPreInversao.assign(rota.begin() + i, rota.begin() + j+1);
            multaRotaPosInversao = multaRotaPreInversao;
            reverse(multaRotaPosInversao.begin(), multaRotaPosInversao.end());

            //Para melhor otimizacao
            //Retiramos a distancia das arestas apagadas e adicionamos as novas
            //Calculamos a nova multa antes da inversao e apos a inversao 
            int novoCusto;
           
            novoCusto = custo - dist[rota[i-1]][rota[i]] - dist[rota[j]][rota[(j+1)%n]] + dist[rota[i-1]][rota[j]] + dist[rota[i]][rota[(j+1)%n]] 
                            - calcDifMulta(i, multaRotaPreInversao , multas) + calcDifMulta(i, multaRotaPosInversao , multas);
            
            vizinho.custo = novoCusto;
            //cout << "O novo custo é " << vizinho.custo << endl;

            vizinhos.push_back(vizinho);

           
            if(otimizar){   // Em casos de grafos grandes podemos adiantar o comeco diminuindo os candidatos
                if(novoCusto<custo){
                    break;
                } 
            }
            
        }
    }
    auto currentTime = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(currentTime - startTime);

    //cout << "Achar os vizinhos gastou " << duration.count() << " segungos" << endl;
    return vizinhos;
}


vector<NewRota> gerarMovimentosVizinhos2(const vector<int>& rota, const int custoRota, const vector<vector<int>> &dist, const vector<vector<int>> &multas,  const vector<vector<int>> &listaTabu) {
    vector<NewRota> vizinhos;
    int n = rota.size();
    //int custo = calcResultado(rota, dist, multas);
    int custo = custoRota;

    //Otimizacao para grafos muito grandes
    bool otimizar=false;
    int contador=0;
    int max = 2000;
    
    auto startTime = std::chrono::steady_clock::now();
    if(n>300){
        otimizar=true;

    }
    

    vector<int> seqi = geraSequenciaAleatoria(1, n-1);
    for (int i = 1; i < n-1; i++) {
        //cout << "Gerei seq" << endl;
        vector<int> seqj = geraSequenciaAleatoria(i+1, n);
        for (int j = i + 1; j < n; j++) {
            if(listaTabu[i][j]!=0) continue;

              if(contador>max){
                return vizinhos;
            }
            if(seqj[j]<seqi[i]) continue;
            //cout << "i: "<< i << "j: " << j << endl;
            //cout << "seqi[i]: " << seqi[i] << " seqj[j]: " << seqj[j] << endl;
           
           
            NewRota vizinho;
            vizinho.rota = rota;
            reverse(vizinho.rota.begin() + seqi[i], vizinho.rota.begin() +seqj[j]+ 1);
            //vizinhos.push_back(vizinho);
            vizinho.i = seqi[i]; vizinho.j = seqj[j];
            //cout <<"Atualizei os vizinho "<< vizinho.i << " " << vizinho.j << endl;
           


            vector<int> multaRotaPreInversao, multaRotaPosInversao;
            //cout << "Antes assign" << endl;
            multaRotaPreInversao.assign(rota.begin() + seqi[i], rota.begin() + seqj[j]+1);
            //cout <<"Gerei rota pre Inv "<< vizinho.i << " " << vizinho.j+1 << endl;
            multaRotaPosInversao = multaRotaPreInversao;
            reverse(multaRotaPosInversao.begin(), multaRotaPosInversao.end());
            //cout <<"fiz REverser "<< endl;

            //Para melhor otimizacao
            //Retiramos a distancia das arestas apagadas e adicionamos as novas
            //Calculamos a nova multa antes da inversao e apos a inversao 
            int novoCusto;
            if(seqi[i]==0){
               novoCusto=calcResultado(rota, dist, multas);
            }else{
                novoCusto = custo - dist[rota[seqi[i]-1]][rota[seqi[i]]] - dist[rota[seqj[j]]][rota[(seqj[j]+1)%n]] + dist[rota[seqi[i]-1]][rota[seqj[j]]] + dist[rota[seqi[i]]][rota[(seqj[j]+1)%n]] 
                                - calcDifMulta(seqi[i], multaRotaPreInversao , multas) + calcDifMulta(seqi[i], multaRotaPosInversao , multas);
            }
            //cout << "calculei novo custo "<< endl;
            vizinho.custo = novoCusto;
            //cout << "O novo custo é " << vizinho.custo << endl;

            vizinhos.push_back(vizinho);

           
            if(otimizar){   // Em casos de grafos grandes podemos adiantar o comeco diminuindo os candidatos
                if(novoCusto<custo){
                    break;
                } 
            }
            contador++;
        }
    }
    auto currentTime = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(currentTime - startTime);

    //cout << "Achar os vizinhos gastou " << duration.count() << " segungos" << endl;
    return vizinhos;
}


void decrementarTabuList(vector<vector<int>> &tabuList){
    int n = tabuList.size();
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            if(tabuList[i][j]!=0){
                tabuList[i][j]-=1;
                tabuList[j][i]-=1;
            } 
                
        }
    }
} 
// Função para executar o algoritmo Tabu Search
vector<int> tabuSearchTDTSP(vector<int> &rota, const vector<vector<int>>& distancias, const vector<vector<int>>& multas, int tamanhoTabu, int maxIteracoes, int maxSegundos) {
    int n = rota.size();
    int maxNoImprove = maxIteracoes/2; //Max iteracoes sem improve
    vector<int> melhorRota = rota;
    int melhorCusto = calcResultado(melhorRota, distancias, multas);
    vector<int> rotaAtual = melhorRota;
    int custoAtual = melhorCusto;

    vector<vector<int>> listaTabu(n, vector<int>(n, 0));
    auto startTime = std::chrono::steady_clock::now();

    for (int iteracao = 0, iteracaoNoImp=0; iteracao < maxIteracoes && iteracaoNoImp<maxNoImprove; iteracao++) {
        auto currentTime = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
        if (duration >= maxSegundos) {  //Maximo x segundos
            break;
        }

      
        
        vector<NewRota> vizinhos = gerarMovimentosVizinhos(rotaAtual, custoAtual,  distancias, multas, listaTabu);
            
        
        
        int melhorVizinhoCusto = numeric_limits<int>::max();
        vector<int> melhorVizinhoRota;
        int i,j;
        for (const auto& vizinho : vizinhos) {
            int vizinhoCusto = vizinho.custo;
            
        

            if (vizinhoCusto < melhorVizinhoCusto) {
                melhorVizinhoCusto = vizinhoCusto;
                melhorVizinhoRota = vizinho.rota;
                i = vizinho.i;
                j = vizinho.j;
                
            }
           
        }
        rotaAtual = melhorVizinhoRota;
        custoAtual = melhorVizinhoCusto;
       //cout << "Atualizei rotas com os valores " << i << " " << j << endl;
        listaTabu[i][j] = listaTabu[j][i] = tamanhoTabu;

        // Atualiza a melhor rota se necessário
        if (melhorVizinhoCusto < melhorCusto) {
            melhorRota = melhorVizinhoRota;
            melhorCusto = melhorVizinhoCusto;
        }else{
            iteracaoNoImp++;
        }
        decrementarTabuList(listaTabu);

    }
    return melhorRota;
}

vector<int> tabuSearch(const vector<vector<int>>& dist, const vector<vector<int>>& multas, int n, int tamanhoTabu, int maxIteracoes, int maxSegundos ){

    vector<int> melhor_rota;
    int melhorCusto = numeric_limits<int>::max();
    auto startTime = std::chrono::steady_clock::now();

    for(int i=0; i<n; i++){
        auto currentTime = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
        if (duration >= maxSegundos) {  //Maximo x segundos
            break;
        }

        // Se o grafo for relativamente pequeno podemos comecar de uma Rota aleatoria
        vector<int> rota = geraRotaAleatoria(dist.size());
        

        rota = tabuSearchTDTSP(rota, dist, multas, tamanhoTabu, maxIteracoes, maxSegundos);
        int newCusto = calcResultado(rota, dist, multas);
        if(newCusto<melhorCusto){
            melhorCusto = newCusto;
            melhor_rota = rota;
        }
    }
    buscaLocalBestImp(dist, melhor_rota, multas);
    return melhor_rota;
}

void grasp(const vector<vector<int>> &dist, vector<int> &rota, const vector<vector<int>> &multas, int maxSegundos, float alp){

    int custoRota = calcResultado(rota,dist, multas);
    vector<int> melhorRota = rota;
    float alpha = alp;
    int maxIterations = 100;
    int n = rota.size();

   
    
    //GRASP
    auto startTime = std::chrono::steady_clock::now();
    while(true){
        auto currentTime = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
        if (duration >= maxSegundos) {
            break;
        }
        
        
        vector<int> novaRota;
        novaRota.push_back(0);
        vector<bool> visitado(rota.size(), false);
        visitado[0]=true;

        for(int i=0; i<n-1; i++){

            int cidade_atual = novaRota[i];

            //auto start1 = std::chrono::steady_clock::now();
            priority_queue<PairWithValue> pq = encontraRLCProximaCidade(dist, visitado,multas, cidade_atual, i+1);
            //auto f2 = std::chrono::steady_clock::now();
            //auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(f2 - start1);
            //cout << "Tempo de execução: " << duration2.count() << " segundos." << "i: "<< i << endl;
            vector<int> elements_value, elements_i;
            // Copiar os elementos da pq para o vector
            while (!pq.empty()) {
                elements_value.push_back(pq.top().value);
                elements_i.push_back(pq.top().key);
                //printf("A cidade %d tem custo %d e ela ",pq.top().key, pq.top().value );
                //cout << visitado[pq.top().key] << endl;
                pq.pop();
            }
            int k = elements_value[0] + alpha*(elements_value[(elements_value.size()-1)]-elements_value[0]);

            //Nossa RLC sera composta de valores menores que k
            int limite_i=n-1;

            for(int j=0; j<elements_value.size(); j++){
                if(elements_value[j]>=k){
                    limite_i = j;
                    //cout << "limite : " << k << endl;
                    break;
                }
            }


            // Escolhemos uma cidade aleatoria
            random_device rd;
            mt19937 rng(rd());
            uniform_int_distribution<int> rdist(0, limite_i);
            int randomNum = rdist(rng);

            cidade_atual = elements_i[randomNum];
            
            novaRota.push_back(cidade_atual);
            visitado[cidade_atual]=true;
           
            
            
            if(n>300)
                buscaLocalBestImpTempo(dist, novaRota, multas, 1);  // Para grafos multo grandes colocamos um limite de tempo na busca local
            else
                buscaLocalBestImp(dist, novaRota, multas);
            
        }
        if(n>300){
            buscaLocalBestImp(dist, novaRota, multas);  //No fim realizamos uma busca local total em n>300
        }


        rota = novaRota;
        //buscaLocalBestImp(dist, rota, multas);
        int custoNovaRota = calcResultado(rota, dist, multas);
        if(custoNovaRota<custoRota){
            custoRota = custoNovaRota;
            melhorRota = rota;
            printf("Achei uma rota melhor!\n");
        }
        //cout << "Gerei uma rota" << endl;
    }
    
    rota = melhorRota;
    buscaLocalBestImp(dist, rota, multas);
}

void buscaLocalBestImp(const vector<vector<int>> &dist, vector<int> &rota, const vector<vector<int>> &multas){  
    bool melhorou = true;
    int n = rota.size();

    int cont=0;
    while(melhorou){
        melhorou = false;
        for(int i=1; i<n-1; i++){   // i=1 garante que a rota[0] nao ira mudar
            for(int j=i+1; j<n; j++){  

                vector<int> multaRotaPreInversao, multaRotaPosInversao;
                multaRotaPreInversao.assign(rota.begin() + i, rota.begin() + j+1);
                multaRotaPosInversao = multaRotaPreInversao;
                reverse(multaRotaPosInversao.begin(), multaRotaPosInversao.end());

                //Para melhor otimizacao
                //Retiramos a distancia das arestas apagadas e adicionamos as novas
                //Calculamos a nova multa antes da inversao e apos a inversao 
                int novoCusto = - dist[rota[i-1]][rota[i]] - dist[rota[j]][rota[(j+1)%n]] + dist[rota[i-1]][rota[j]] + dist[rota[i]][rota[(j+1)%n]] 
                                - calcDifMulta(i, multaRotaPreInversao , multas) + calcDifMulta(i, multaRotaPosInversao , multas);
                
                //Caso novoCusto seja negativo quer dizer que haverá um  ganho, entao modificamos a rota
                if(novoCusto < 0){
                    melhorou = true;
                    reverse(rota.begin() + i, rota.begin() + j + 1 ); //Nos "apagamos" os arcos (i,i+1) e o arco(j+1,j+2)
                }

            }
        }
    }

    return;

}

void buscaLocalBestImp2(const vector<vector<int>> &dist, vector<int> &rota, const vector<vector<int>> &multas){  //Versao 2-opt Antiga (má otimizada), utilize apenas para conferir resultados/testes
    int custo = calcResultado(rota, dist, multas);
    bool melhorou = true;
    int n = rota.size();

    int cont=0;
    while(melhorou){
        melhorou = false;
        for(int i=1; i<n-1; i++){   // i=1 garante que a rota[0] nao ira mudar
            for(int j=i+1; j<n; j++){  

                vector<int> rotaAux = rota;
                reverse(rotaAux.begin() + i, rotaAux.begin() + j + 1); //Nos "apagamos" os arcos (i,i+1) e o arco(j+1,j+2)
                int novoCusto = calcResultado(rotaAux, dist, multas);  //Custo da rotaAux

                if(novoCusto < custo){  //Caso o novoCusto seja melhor, rota = rotaAux
                    custo = novoCusto;
                    melhorou = true;
                    rota = rotaAux;
                }

            }
        }
    }

    return;

}

void buscaLocalBestImpTempo(const vector<vector<int>> &dist, vector<int> &rota, const vector<vector<int>> &multas, int maxSegundos){  
    bool melhorou = true;
    int n = rota.size();

    int cont=0;
    auto startTime = std::chrono::steady_clock::now();
 
        
        
    while(melhorou){
        melhorou = false;
        for(int i=1; i<n-1; i++){   // i=1 garante que a rota[0] nao ira mudar
            for(int j=i+1; j<n; j++){  
                auto currentTime = std::chrono::steady_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
                if (duration >= maxSegundos) {
                    return;
                }

                vector<int> multaRotaPreInversao, multaRotaPosInversao;
                multaRotaPreInversao.assign(rota.begin() + i, rota.begin() + j+1);
                multaRotaPosInversao = multaRotaPreInversao;
                reverse(multaRotaPosInversao.begin(), multaRotaPosInversao.end());

                //Para melhor otimizacao
                //Retiramos a distancia das arestas apagadas e adicionamos as novas
                //Calculamos a nova multa antes da inversao e apos a inversao 
                int novoCusto = - dist[rota[i-1]][rota[i]] - dist[rota[j]][rota[(j+1)%n]] + dist[rota[i-1]][rota[j]] + dist[rota[i]][rota[(j+1)%n]] 
                                - calcDifMulta(i, multaRotaPreInversao , multas) + calcDifMulta(i, multaRotaPosInversao , multas);
                
                //Caso novoCusto seja negativo quer dizer que haverá um  ganho, entao modificamos a rota
                if(novoCusto < 0){
                    melhorou = true;
                    reverse(rota.begin() + i, rota.begin() + j + 1 ); //Nos "apagamos" os arcos (i,i+1) e o arco(j+1,j+2)
                }

            }
        }
    }
    
    return;

}


void leDados(vector<pair<float,float>> &coord, vector<vector<int>> &matriz){  
   
    int qtd_linhas = coord.size();
    for(int i=0; i<qtd_linhas; i++){ 
        int aux;
        cin>>aux; // Ignorara numero da cidade pois ja temos essa informacao

        cin >> coord[i].first;
        cin >> coord[i].second;
    }
    
    //Lera Multas
    



    int qtd_m;
    cin >> qtd_m;   //Comentar caso nao tenha arquiivo de multas 0
    for(int i=0; i<qtd_linhas; i++){
        for(int j=0; j<qtd_linhas; j++){
            //matriz[i][j]=0;  //Descomentar caso nao tenho arquivo de multas 0 e comentar cin abaixo
            cin >> matriz[i][j];
        }
    }
    
    return;
}

void leDadosZero(vector<pair<float,float>> &coord, vector<vector<int>> &matriz){  //Utilize essa funcao caso instancia nao tenha .txt para multas 0
   
    int qtd_linhas = coord.size();
    for(int i=0; i<qtd_linhas; i++){ 
        int aux;
        cin>>aux; // Ignorara numero da cidade pois ja temos essa informacao

        cin >> coord[i].first;
        cin >> coord[i].second;
    }
    
    //Lera Multas
    

    int qtd_m;
    //cin >> qtd_m;   //Comentar caso nao tenha arquivo de multas 0
    for(int i=0; i<qtd_linhas; i++){
        for(int j=0; j<qtd_linhas; j++){
            matriz[i][j]=0;  //Descomentar caso nao tenho arquivo de multas 0 e comentar cin abaixo
            //cin >> matriz[i][j];
        }
    }
    
    return;
}

int main(){

    int qtd_linhas;
    cin >> qtd_linhas;
    vector<pair<float,float>> coord(qtd_linhas);    //Vetor que guardara coordenadas das cidades
    vector<vector<int>> multas(qtd_linhas, vector<int>(qtd_linhas));    //Guarda as multas da cidades

    leDados(coord, multas); //Le os dados
    //leDadosZero(coord, multas); //Caso instancia nao tenho .txt para multas 0
    
    auto inicio = std::chrono::high_resolution_clock::now(); // Obtém o tempo atual

    vector<vector<int>> dist = geraMatrizDist(coord);
    //vector<int> rota = geraRotaSequencial(coord.size());
    vector<int> rota = geraRotaAleatoria(coord.size());
    //vector<int> rota = geraRotaGuloso(dist, coord.size());
    //buscaLocalBestImp(dist, rota, multas);

    printRota(rota);
    cout << endl;
    //buscaLocalBestImp(dist, rota, multas);

    rota = tabuSearch( dist, multas,30,  (1/2)*rota.size(), 2000, 5*60);
    //grasp(dist, rota, multas, 2*60, 0.4);
 
    printRota(rota);

    auto fim = std::chrono::high_resolution_clock::now(); // Obtém o tempo após a execução
    auto duracao = std::chrono::duration_cast<std::chrono::seconds>(fim - inicio); // Calcula a diferença em microssegundos

    cout << "Tempo de execução: " << duracao.count() << " segundos." << endl;
    
    cout << "O resultado e " << calcResultado(rota, dist, multas) << endl;

    
    return 0;
}