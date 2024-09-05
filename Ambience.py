import numpy as np
from CoolProp.HumidAirProp import HAPropsSI
import math


# declaração de classes
class Ambiente:
    def __init__(self, nome, t, phi, n_pessoas, area):
        self.nome = nome
        self.t = t  # temperatura do ambiente
        self.phi = phi  # umidade relativa
        self.n_pessoas = n_pessoas
        self.area = area
        self.Q_sens = 0  # calculado depois
        self.Q_lat = 0
        self.alfa = 0.5
        self.V_pessoas = n_pessoas * Vpp
        self.V_renov = self.area * pdir * 10 / 3600
        self.V = max(self.V_pessoas, self.V_renov)
        self.q_cont = 0


class ParedeExt:
    def __init__(self, pontos_s, orient, amb_in, material, cor):
        self.pontos_s = pontos_s  # lista de pontos, indice ou coordenadas, não sei ainda, decido depois
        self.orient = orient  # norte, sul, leste, oeste
        self.amb_in = amb_in  # ambiente interno a esta parede
        self.material = material
        self.cor = cor
        self.area = None  # calculado depois por uma função
        self.u = U[material]
        self.alfa = Alfa[cor]
        self.I = I[orient]


class ParedeInt:
    def __init__(self, pontos_s, material, amb_in, amb_ex):
        self.pontos_s = pontos_s
        self.material = material
        self.amb_in = amb_in
        self.amb_ex = amb_ex
        self.area = None
        self.u = U[material]


class Transparente:
    def __init__(self, parede, A, n, material, protecao):
        self.parede = parede  # em qual parede a janela fica
        self.A = A  # area
        self.n = n  # número de objetos iguais a este
        self.material = material
        self.protecao = protecao
        self.fs1 = Fs1[material]
        self.fs2 = Fs2[protecao]


class FonteCalor:
    def __init__(self, ambiente, n, w, eta, ft, fu):
        self.ambiente = ambiente
        self.n = n
        self.w = w
        self.eta = eta
        self.ft = ft
        self.fu = fu


class FonteContaminante:
    def __init__(self, ambiente, n, q):
        self.ambiente = ambiente
        self.n = n
        self.q = q


class Iluminacao:
    def __init__(self, ambiente, n, q, fu, fc, ft):
        self.ambiente = ambiente
        self.n = n
        self.q = q
        self.fu = fu
        self.fc = fc
        self.ft = ft


class Dutos_Segmentos:
    def __init__(self, inicio, fim, bool_primeiro):
        self.inicio = inicio
        self.fim = fim
        self.bool_primeiro = bool_primeiro
        self.L = math.sqrt(
            (fim.coords[0] - inicio.coords[0]) ** 2 + (fim.coords[1] - inicio.coords[1]) ** 2)  # extensão da seção
        self.u = None  # velocidade de fluxo
        self.delta_Q = self.inicio.Q  # vazão volumétrica para a boca
        self.Q = 0
        self.A_st = None  # self.Q/self.u  # área da seção transversal
        self.h = None
        self.l = None  # largura da seção transversal
        self.delta_P = None  # 0.001026 * self.L * self.u ** 2.5 / self.Q ** 0.61


class Dutos_Elementos:
    def __init__(self, F, coords):
        self.F = F
        self.coords = coords
        self.u = None
        self.delta_P = None  # F * rho_ar * self.u ** 2 / (2 * 9.81)
        self.Q = 0  # elementos de dutos nâo consomem fluxo volumétrico, ao contrário de bocas de insuflamento


class Grelhas:
    def __init__(self, div, Q, T, coords):
        self.div = div
        self.Q = Q
        self.T = T
        self.coords = coords
        self.a = aGrelhas[div]
        self.K = KGrelhas[div]
        self.omega = 1 / self.a * (self.K * Q / T) ** 2
        self.cf = Q / self.omega
        self.c = self.cf / self.a
        self.delta_P = 0


# VARIÁVEIS GLOBAIS

pdir = 3.0  # pé direito
rho_ar = 1.2924  # densidade do ar, kg/m**3
u_telhado = 2.04  # laje concreto 10cm + fibrocimento. O trabalho cita só o fibrocimento, mas o professor disse para
# usar fibrocimento com concreto porque não faria sentido ter só o fibrocimento
qsp = 83.84  # calor sensível por pessoa, dados do exemplo
qlp = 32.56  # calor latente por pessoa, dados do exemplo
Vpp = 13 / 3600  # fluxo volumétrico de ar por pessoa, dados do exemplo
A_total = 20.0 * 30.0  # área do prédio
A_tridimensionais = 5.0 * 15.0  # área da sala das tridimensionais
A_banheiro = 3.0 * 10.0  # área do banheiro
A_principal = A_total - A_banheiro - A_tridimensionais  # área da sala principal
R_t = 0.04  # resistência térmica superficial

# dicionários para extrair propriedades dos materiais e parâmetros
Alfa = {
    "escura": 0.9,
    "média": 0.7,
    "clara": 0.5
}
U = {
    "tijolos_15cm": 2.88 * 1.1627, #kcal/h -> W,
    "drywall": 4.26,
    "vidro": 5.18 * 1.1627,
    "fumê": 5.18 * 1.1627
}
I = {
    "Sul": 113.5,
    "Leste": 551.8,
    "Norte": 56.7,
    "Oeste": 510.8
}
Ambientes = {
    "externo": 0,
    "principal": 1,
    "tridimensionais": 2,
    "banheiro": 3
}
Fs1 = {
    "vidro": 0.87,  # transparente simples 3mm
    "fumê": 0.6  # fumê 6mm
}
Fs2 = {
    "persiana": 0.6,
    "sem": 1
}
KGrelhas = {
    30: 9.5,
    45.2: 7.6,
    60: 6.2,
    90: 5
}
aGrelhas = {
    30: 0.68,
    45.2: 0.62,
    60: 0.62,
    90: 0.58
}

# U de tijolos 15cm: 2.88
# U de drywall 12.5mm: 4.26 / fonte: https://www.otm.sg/u-value-calculator, espessura: busca por drywall comercial

# listas de cada instância de cada classe
ambientes = []
paredes_ext = []
paredes_int = []
transparentes = []
fontes_calor = []
pessoas = []
fontes_contaminante = []
iluminacao = []
dutos = [[],[],[]]
grelhas = []
elementos = []

# declaração dos ambientes. ASSUMIDO: Calculando para a situação crítica de cada ambiente, então considera-se que
# todos os ambientes estão na lotação máxima simultaneamente. Isso soma mais pessoas do que de fato trabalham no
# laboratório, mas garante que os parâmetros de condicionamento são atingidos em todas as situações.
externo = Ambiente("externo", 40, 80, 0, 999999999)  # área do ambiente externo não importa
principal = Ambiente("principal", 23, 40, 10, A_principal)
tridimensionais = Ambiente("tridimensionais", 20, 30, 10, A_tridimensionais)
banheiro = Ambiente("banheiro", 25, 50, 4, A_banheiro)

ambientes.extend([externo, principal, tridimensionais, banheiro])

# longo bloco de declaração de propriedades das paredes

pontos_paredes_coords = [[0, 0], [0, 20], [27, 20], [30, 20], [30, 10], [30, 5], [30, 0], [15, 0], [15, 5], [27, 10]]
pontos_paredes_externas_pontas = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 0]]
pontos_paredes_internas_pontas = [[7, 8], [8, 5], [9, 4], [9, 2]]

# todos os delta_T equivalentes de radiação serão calculados

# passa para cada objeto Parede o objeto Ambiente correspondente. Acho que dá pra tirar as temperaturas assim?

# em cada lista, associa o ambiente interno das paredes internas com os índices especificados
parede_ext_ambiente_int = []
parede_ext_ambiente_int_principal = [0, 1, 4, 7]
parede_ext_ambiente_int_banheiro = [2, 3]
parede_ext_ambiente_int_tridimensionais = [5, 6]

# constroi a lista de ambientes internos de cada parede externa
for i in range(len(pontos_paredes_externas_pontas)):
    flag_interior = 0
    for j in range(len(parede_ext_ambiente_int_principal)):
        if i == parede_ext_ambiente_int_principal[j]:
            parede_ext_ambiente_int.append(principal)
            flag_interior = 1
            break
    for j in range(len(parede_ext_ambiente_int_banheiro)):
        if i == parede_ext_ambiente_int_banheiro[j]:
            parede_ext_ambiente_int.append(banheiro)
            flag_interior = 1
            break
    for j in range(len(parede_ext_ambiente_int_tridimensionais)):
        if i == parede_ext_ambiente_int_tridimensionais[j]:
            parede_ext_ambiente_int.append(tridimensionais)
            flag_interior = 1
            break
    if flag_interior == 0:
        print("Índice da parede:", i)
        raise Exception("Parede externa sem ambiente interior definido")

parede_int_ambiente_ext = []
for i in range(len(pontos_paredes_internas_pontas)):
    parede_int_ambiente_ext.append(principal)

parede_int_ambiente_int = []  # np.zeros(len(pontos_paredes_internas_pontas))
parede_int_ambiente_int_tridimensionais = [0, 1]
parede_int_ambiente_int_banheiro = [2, 3]
for i in range(len(pontos_paredes_internas_pontas)):
    flag_interior = 0
    for j in range(len(parede_int_ambiente_int_tridimensionais)):
        if i == parede_int_ambiente_int_tridimensionais[j]:
            parede_int_ambiente_int.append(tridimensionais)
            flag_interior = 1
    for j in range(len(parede_int_ambiente_int_banheiro)):
        if i == parede_int_ambiente_int_banheiro[j]:
            parede_int_ambiente_int.append(banheiro)
            flag_interior = 1
    if flag_interior == 0:
        print("Índice da parede:", i)
        raise Exception("Parede interna sem ambiente interior definido")

# parede_ext_orientacao = np.zeros(len(pontos_paredes_externas_pontas))
parede_ext_orientacao = []
parede_ext_orientacao_norte = [1, 2]
parede_ext_orientacao_sul = [6, 7]
parede_ext_orientacao_leste = [3, 4, 5]
parede_ext_orientacao_oeste = [0]

for i in range(len(pontos_paredes_externas_pontas)):
    flag_interior = 0
    for j in range(len(parede_ext_orientacao_norte)):
        if i == parede_ext_orientacao_norte[j]:
            parede_ext_orientacao.append("Norte")
            flag_interior = 1
    for j in range(len(parede_ext_orientacao_sul)):
        if i == parede_ext_orientacao_sul[j]:
            parede_ext_orientacao.append("Sul")
            flag_interior = 1
    for j in range(len(parede_ext_orientacao_leste)):
        if i == parede_ext_orientacao_leste[j]:
            parede_ext_orientacao.append("Leste")
            flag_interior = 1
    for j in range(len(parede_ext_orientacao_oeste)):
        if i == parede_ext_orientacao_oeste[j]:
            parede_ext_orientacao.append("Oeste")
            flag_interior = 1
    if flag_interior == 0:
        print("Uma parede externa está sem orientacao definida. Índice da parede:", i)
        raise Exception("Parede externa sem orientacao definida")

material_paredes_int = []
material_paredes_ext = []
cor_paredes_ext = []

for i in range(len(pontos_paredes_externas_pontas)):
    material_paredes_ext.append("tijolos_15cm")
    cor_paredes_ext.append("clara")
for i in range(len(pontos_paredes_internas_pontas)):
    material_paredes_int.append("drywall")

# construção de paredes
for i in range(len(pontos_paredes_externas_pontas)):
    paredes_ext.append(ParedeExt(pontos_paredes_externas_pontas[i],
                                 parede_ext_orientacao[i],
                                 parede_ext_ambiente_int[i],
                                 material_paredes_ext[i],
                                 cor_paredes_ext[i], ))
for i in range(len(pontos_paredes_internas_pontas)):
    paredes_int.append(ParedeInt(pontos_paredes_internas_pontas[i],
                                 material_paredes_int[i],
                                 parede_int_ambiente_int[i],
                                 parede_int_ambiente_ext[i]))
# com isso, todas as propriedades das paredes estão definidas. os outros objetos vou declarar instância por instância.

janela_banheiro_norte = Transparente(2, 1 * 0.5, 1, "vidro", "persiana")
janela_banheiro_leste = Transparente(3, 1 * 0.5, 1, "vidro", "persiana")
janela_3d_sul = Transparente(6, 2 * 1, 2, "vidro", "persiana")
janela_principal_sul = Transparente(7, 2 * 1, 3, "vidro", "persiana")
porta = Transparente(0, 3 * 3, 1, "fumê", "sem")
transparentes.extend([janela_banheiro_norte, janela_banheiro_leste, janela_3d_sul, janela_principal_sul, porta])

fonte_tridimensional = FonteCalor(tridimensionais, 4, 300, 1, 1, 1)
fontes_calor.append(fonte_tridimensional)

benzina = FonteContaminante(principal, 1, 0.174)  # q é dado em m**3/h e é tirado do exemplo no moodle
fontes_contaminante.append(benzina)

iluminacao_principal = Iluminacao(principal, 200, 65 + 16, 1, 1.2, 1)
iluminacao_tridimensionais = Iluminacao(tridimensionais, 30, 65 + 16, 1, 1.2, 1)
iluminacao_banheiro = Iluminacao(banheiro, 10, 65 + 16, 1, 1.2, 1)
iluminacao.extend([iluminacao_principal, iluminacao_tridimensionais, iluminacao_banheiro])

for i in range(len(pontos_paredes_externas_pontas)):
    coord_x_1 = pontos_paredes_coords[paredes_ext[i].pontos_s[0]][0]
    coord_x_2 = pontos_paredes_coords[paredes_ext[i].pontos_s[1]][0]
    coord_y_1 = pontos_paredes_coords[paredes_ext[i].pontos_s[0]][1]
    coord_y_2 = pontos_paredes_coords[paredes_ext[i].pontos_s[1]][1]
    area_prelim = np.sqrt((coord_x_2 - coord_x_1) ** 2 + (coord_y_2 - coord_y_1) ** 2) * pdir
    area_jan = 0
    for j in range(len(transparentes)):
        if transparentes[j].parede == i:
            area_jan += transparentes[j].A * transparentes[j].n
    paredes_ext[i].area = area_prelim - area_jan

for i in range(len(pontos_paredes_internas_pontas)):
    coord_x_1 = pontos_paredes_coords[paredes_int[i].pontos_s[0]][0]
    coord_x_2 = pontos_paredes_coords[paredes_int[i].pontos_s[1]][0]
    coord_y_1 = pontos_paredes_coords[paredes_int[i].pontos_s[0]][1]
    coord_y_2 = pontos_paredes_coords[paredes_int[i].pontos_s[1]][1]
    paredes_int[i].area = np.sqrt((coord_x_2 - coord_x_1) ** 2 + (coord_y_2 - coord_y_1) ** 2) * pdir

# ======== BLOCO DE CÁLCULOS ==========

# carga térmica para cada ambiente

# [ I M P O R T A N T E ]

# ADMITE-SE que o calor que entra em um ambiente por penetração vindo de uma parede interna SUBTRAI da carga térmica do
# ambiente de maior temperatura

for i in range(1,
               len(ambientes)):  # 0 é o ambiente externo; não faz sentido calcular a carga térmica para o ambiente externo, por isso parte de 1

    # penetração

    print(ambientes[i].nome, ":")

    Qp_ext = 0
    for j in range(len(paredes_ext)):
        if paredes_ext[j].amb_in == ambientes[i]:
            Qp_ext += paredes_ext[j].area * paredes_ext[j].u * (
                        ambientes[0].t - paredes_ext[j].amb_in.t)  # penetração vinda de paredes externas
    ambientes[i].Q_sens += Qp_ext
    print("Qp_ext = ", Qp_ext)

    Qp_int = 0
    for j in range(len(paredes_int)):
        if paredes_int[j].amb_in == ambientes[i]:
            Qp_int += paredes_int[j].area * paredes_int[j].u * (
                        paredes_int[j].amb_ex.t - paredes_int[j].amb_in.t)  # penetração vinda de paredes internas
    ambientes[i].Q_sens += Qp_int  # Q_int não necessariamente é positivo
    print("Qp_int = ", Qp_int)

    # O for acima calcula o calor de penetração advindo das paredes internas CUJOS AMBIENTES INTERNOS sejam ambientes[i].
    # abaixo, será feito o mesmo for, mas para o calor advindo de paredes internas cujos ambientes EXTERNOS sejam ambientes[i]
    # para uma mesma parede, isso dá o mesmo valor, porém com o sinal trocado. Conforme suposto, o calor de
    # penetração que adentra um ambiente vindo de uma parede interna deve sair de outro ambiente

    Qp_int = 0
    for j in range(len(paredes_int)):
        if paredes_int[j].amb_ex == ambientes[i]:
            Qp_int += paredes_int[j].area * paredes_int[j].u * (
                        paredes_int[j].amb_in.t - paredes_int[j].amb_ex.t)  # penetração vinda de paredes internas
    ambientes[i].Q_sens += Qp_int
    print("Qp_int = ", Qp_int)

    # penetração pelo telhado

    ambientes[i].Q_sens += ambientes[i].area * u_telhado * (ambientes[0].t - ambientes[i].t)
    print("Qp_telhado = ", ambientes[i].area * u_telhado * (ambientes[0].t - ambientes[i].t))

    # insolação opaca

    Qi_op = 0
    for j in range(len(paredes_ext)):  # só ocorre em paredes externas
        if paredes_ext[j].amb_in == ambientes[i]:
            delta_t_eq = paredes_ext[j].alfa * R_t * paredes_ext[j].I
            Qi_op += paredes_ext[j].area * paredes_ext[j].u * delta_t_eq
    ambientes[i].Q_sens += Qi_op
    print("Qi_op = ", Qi_op)

    # insolação transparente

    Qi_tr = 0
    for j in range(len(paredes_ext)):  # para cada parede externa
        if paredes_ext[j].amb_in == ambientes[i]:  # se o ambiente interno a essa parede for ambiente[i]
            for k in range(len(transparentes)):  # para cada elemento transparente
                if transparentes[k].parede == j:  # se esse elemento pertence a parede_ext[j]
                    Qi_tr += transparentes[k].A * transparentes[k].fs1 * transparentes[k].fs2 * paredes_ext[
                        j].I  # calcula Qi_tr
    ambientes[i].Q_sens += Qi_tr
    print("Qi_tr = ", Qi_tr)

    # insolação pelo telhado

    delta_t_eq_telh = ambientes[i].alfa * R_t * 857.8
    ambientes[i].Q_sens += ambientes[i].area * u_telhado * delta_t_eq_telh
    print("Qi_telhado = ", ambientes[i].area * u_telhado * delta_t_eq_telh)
    # iluminação

    Q_il = 0
    for j in range(len(iluminacao)):
        if iluminacao[j].ambiente == ambientes[i]:
            Q_il += iluminacao[j].n * iluminacao[j].q * iluminacao[j].fu * iluminacao[j].fc * iluminacao[j].ft
    ambientes[i].Q_sens += Q_il
    print("Q_il = ", Q_il)

    # pessoas

    ambientes[i].Q_sens += ambientes[i].n_pessoas * qsp
    print("Q_sens_pess = ", ambientes[i].n_pessoas * qsp)
    ambientes[i].Q_lat += ambientes[i].n_pessoas * qlp
    print("Q_lat_pess = ", ambientes[i].n_pessoas * qlp)

    # equipamentos

    Q_eq = 0
    for j in range(len(fontes_calor)):
        if fontes_calor[j].ambiente == i:
            Q_eq = fontes_calor[j].w / fontes_calor[j].eta * fontes_calor[j].ft * fontes_calor[j].fu
    ambientes[i].Q_sens += Q_eq
    print("Q_eq = ", Q_eq)

    # ar externo

    Cp = HAPropsSI("cp_ha", "T", ambientes[i].t + 273.15, "R", ambientes[i].phi / 100, "P", 101325)
    print(Cp)
    omega = HAPropsSI("Omega", "T", ambientes[i].t + 273.15, "R", ambientes[i].phi / 100, "P", 101325)
    omega_ext = HAPropsSI("Omega", "T", ambientes[0].t + 273.15, "R", ambientes[0].phi / 100, "P", 101325)
    hg = HAPropsSI("H", "T", ambientes[i].t + 273.15, "R", ambientes[i].phi / 100, "P", 101325)
    hg_ext = HAPropsSI("H", "T", ambientes[i].t + 273.15, "R", ambientes[i].phi / 100, "P", 101325)
    ambientes[i].Q_sens += ambientes[i].V * 1.225 * Cp * (ambientes[0].t - ambientes[i].t)
    print("Q_sens_ar = ", ambientes[i].V * 1.225 * Cp * (ambientes[0].t - ambientes[i].t))
    ambientes[i].Q_lat += ambientes[i].V * 1.225 * (omega_ext * hg_ext - omega * hg)
    print("Q_lat_ar = ", ambientes[i].V * 1.225 * (omega_ext * hg_ext - omega * hg))

    print(ambientes[i].nome, ": resumo")
    print("Q_sens total = ", ambientes[i].Q_sens)
    print("Q_lat total = ", ambientes[i].Q_lat)
    print("V = ", ambientes[i].V)
    print("=======")


# ===================== CONDICIONAMENTO DE AR ================

# primeiro ambiente: principal
# arbitrando a temperatura de condicionamento como a mínima admissível
# T_cond = ambientes[1].t - 13
# V_cond = ambientes[1].Q

# caso crítico: sala das tridimensionais

# um_abs = HAPropsSI("RelHum", "T", ambientes[2].t + 273.15, "R", ambientes[2].phi / 100, "P", 101325)
# h = HAPropsSI("H", "T", ambientes[2].t + 273.15, "R", ambientes[2].phi / 100, "P", 101325)




# ==================== PERDA DE CARGA =====================

grelhas.extend([Grelhas(60, ambientes[1].V / 7 * 1.2, 10, [10.7734, 10]),  # 0; ambientes[1] = principal
                Grelhas(60, ambientes[1].V / 7 * 1.2, 10, [21.2265, 10]),  # 1
                Grelhas(60, ambientes[1].V / 7 * 0.8, 7.5, [7.5, 15.6699]),  # 2
                Grelhas(90, ambientes[1].V / 7 * 1.2, 7.5, [7.5, 7.5]),  # 3
                Grelhas(90, ambientes[1].V / 7 * 1.8, 10, [10, 10]),  # 4
                Grelhas(45.2, ambientes[1].V / 7 * 0.2, 5, [17.5, 10]),  # 5
                Grelhas(90, ambientes[1].V / 7 * 0.6, 5, [25, 10]),  # 6
                Grelhas(60, ambientes[2].V / 3, 5, [17.8868, 0]),  # 7; ambientes[2] = tridimensionais
                Grelhas(60, ambientes[2].V / 3, 5, [25.5018, 0]),  # 8
                Grelhas(60, ambientes[2].V / 3, 5, [27.1132, 0]),  # 9
                Grelhas(90, ambientes[3].V / 2, 3, [27, 13]),  # 10; ambientes[3] = banheiro
                Grelhas(90, ambientes[3].V / 2, 3, [27, 17])])  # 11

for i in range(len(grelhas)):
    if grelhas[i].c > 2.5 and grelhas[i].c < 3.8:
        print("i", i, ": c = ", grelhas[i].c, "m/s")
    else:
        print("Índice:",i,"; velocidade:", grelhas[i].c, "m/s")
        raise Exception("Grelha com velocidade de saída fora dos parâmetros")
    print("omega = ", grelhas[i].omega)
    grelhas[i].delta_P = 1 * (grelhas[i].c) ** 2 * rho_ar / (2 * 9.81)  # considerando lambda = 1

elementos.extend([Dutos_Elementos(0.92, [27, 10]),
                  # para cotovelos: todos com ângulo de 90º, l/h = 3 e r = 0, portanto F=0.92
                  Dutos_Elementos(1, [7.5, 10]),  # para ramal: ângulo de 90º, F=1
                  Dutos_Elementos(0.92, [7.5, 0]),
                  Dutos_Elementos(0,[7.5,20])]) # elemento falso, apenas para contabilizar a perda de carga até o CAE

# poderia implementar uma maneira de
# identificar os ramos ao invès de fazer
# três listas hard-coded, mas infelizmente
# meu tempo e minha sanidade estão no fim

dutos[0].extend([Dutos_Segmentos(grelhas[11], grelhas[10], True),
              Dutos_Segmentos(grelhas[10], elementos[0], False),
              Dutos_Segmentos(elementos[0], grelhas[6], False),
              Dutos_Segmentos(grelhas[6], grelhas[1], False),
              Dutos_Segmentos(grelhas[1], grelhas[5], False),
              Dutos_Segmentos(grelhas[5], grelhas[0], False),
              Dutos_Segmentos(grelhas[0], grelhas[4], False),
              Dutos_Segmentos(grelhas[4], elementos[1], False)])

dutos[1].extend([Dutos_Segmentos(grelhas[9], grelhas[8], True),
                  Dutos_Segmentos(grelhas[8], grelhas[7], False),
                  Dutos_Segmentos(grelhas[7], elementos[2], False),
                  Dutos_Segmentos(elementos[2], grelhas[3], False),
                  Dutos_Segmentos(grelhas[3],elementos[1], False)])

soma_Q = 0 #soma os fluxos na divergência dos ramais

for i in range(2):
    for j in range(len(dutos[i])):
        if dutos[i][j].bool_primeiro == True:
            dutos[i][j].u = 3.9
        dutos[i][j].Q += dutos[i][j].delta_Q
        dutos[i][j].A_st = dutos[i][j].Q / dutos[i][j].u
        dutos[i][j].h = math.sqrt(dutos[i][j].A_st / 3)  # dividido por 3 para que l/h = 3
        dutos[i][j].l = dutos[i][j].A_st / dutos[i][j].h
        print("segmento entre",dutos[i][j].inicio.coords,"e",dutos[i][j].fim.coords,":")
        print("h =", dutos[i][j].h,"m")
        print("l =", dutos[i][j].l, "m")
        dutos[i][j].delta_P = 0.001026 * dutos[i][j].L * dutos[i][j].u ** 2.5 / dutos[i][j].Q ** 0.61
        if j < len(dutos[i]) - 1:
            dutos[i][j + 1].u = math.sqrt(dutos[i][j].u ** 2 + dutos[i][j].delta_P * 2 * 9.81 / (0.75 * rho_ar))
            dutos[i][j + 1].Q = dutos[i][j].Q
    soma_Q += dutos[i][j].Q


# específico para a última seção, diretamente à jusante do CAE

dutos[2].extend([Dutos_Segmentos(elementos[1], grelhas[2], True),
                  Dutos_Segmentos(grelhas[2], elementos[3], False)])

dutos[2][0].Q = soma_Q # o fluxo no início desta seção é a soma dos ramais à jusante

for j in range(len(dutos[2])):
    if dutos[2][j].bool_primeiro == True:
        dutos[2][j].u = 3.9
    dutos[2][j].Q += dutos[2][j].delta_Q
    dutos[2][j].A_st = dutos[2][j].Q / dutos[2][j].u
    dutos[2][j].h = math.sqrt(dutos[2][j].A_st / 3)  # dividido por 3 para que l/h = 3
    dutos[2][j].l = dutos[2][j].A_st / dutos[2][j].h
    dutos[2][j].delta_P = 0.001026 * dutos[2][j].L * dutos[2][j].u ** 2.5 / dutos[2][j].Q ** 0.61
    if j < len(dutos[2]) - 1:
        dutos[2][j + 1].u = math.sqrt(dutos[2][j].u ** 2 + dutos[2][j].delta_P * 2 * 9.81 / (0.75 * rho_ar))
        dutos[2][j + 1].Q = dutos[2][j].Q

# PERDA DE CARGA NO CAE

V_CAE = 0
for i in range(1,len(ambientes)):
    V_CAE += ambientes[i].V
u_CAE = 4.5 # m/s, velocidade de entrada, dos slides
A_CAE = V_CAE/u_CAE
a_CAE = 0.75
lambda_CAE = 2
c_CAE = u_CAE/a_CAE
delta_P_CAE = lambda_CAE * c_CAE ** 2 * rho_ar / (2 * 9.81)

# PERDA DE CARGA TOTAL

delta_P_total = 0
for i in range(len(grelhas)):
    delta_P_total += grelhas[i].delta_P
for i in range(len(dutos)):
    for j in range(len(dutos[i])):
        delta_P_total += dutos[i][j].delta_P
delta_P_total += delta_P_CAE

print("===========")
print("DELTA P TOTAL:", round(delta_P_total,3), "mmH2O,", round(delta_P_total * 997 * 9.81 * 10 ** -3, 3) , "Pa")
print("POTÊNCIA DE VENTILAÇÃO:", round(delta_P_total * 997 * 9.81 * 10 ** -3 * V_CAE, 3), "W")