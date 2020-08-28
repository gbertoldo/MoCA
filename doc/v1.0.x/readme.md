MoCA v1.0.0

MoCA - Method of Characteristics for Axisymmetric flows

Breve descrição
O código MoCA resolve o escoamento compressível de gases através de tubeiras de simetria axial utilizando o método das características [??]. Com relação ao perfil da tubeira, há duas possibilidades: pode ser prescrito pelo usuário ou pode ser determinada pelo método de Rao[??] para que maximize o coeficiente de empuxo. Os resultados produzidos incluem o coeficiente de empuxo, a rede de características com a distribuição do número de Mach, do ângulo do vetor velocidade com relação à direção axial.

Como compilar o código fonte

Pré-requisitos para compilação:
- compilador para C++
- biblioteca científica GNU (GSL - GNU Scientific Library)


Como usar o programa

Basta executar o seguinte comando
MoCA <config.json>
onde config.json é o nome do arquivo de configuração. Este arquivo deve estar no formato json??.

Parâmetros de configuração


Linha inicial prescrita através de um arquivo de texto simples ou gerada através da solução de Kliegel-Levine.

Perfil do convergente:
Perfil do divergente:

