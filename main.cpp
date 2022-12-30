#include <iostream>
#include <random>
#include <vector>
#include <cmath>

#include <SFML/Graphics.hpp>

#define MARGIN 20.0f
#define WINDOW_HEIGHT 600
#define WINDOW_WIDTH 800

#define VIEW_HEIGHT (float(WINDOW_HEIGHT)-2*MARGIN)
#define VIEW_WIDTH (float(WINDOW_WIDTH)-2*MARGIN)

double LY = 6.0;
double LX = 3.0;

#define V_SCALE
#ifdef V_SCALE
float factor = VIEW_HEIGHT / float(LY);
#else
float factor = VIEW_WIDTH / float(LX);
#endif

double H0 = 6.0;
#define NY 25
#define NX 15
double DY = LY / double(NY);
double DX = LX / double(NX);

double W = 3.0;
double H = 2.0;

// PHYSICAL CONSTANTS
// Gravitational acceleration [m/s^2]
double g = 9.81;
// Particle diameter [m]
double D = 0.15;
// Spring stiffness [N/m]
double kc = 10000.0;
// Coefficient of restitution
double en = 0.5;
// Particle mass [kg]
double M = 0.1;
// Damping coefficient
double eta = -2.0*sqrt(M*kc)*log(en)/sqrt(log(en)*log(en) + M_PI*M_PI);

// Maps physical coordinates to sfml coordinates
sf::Vector2f phys2sf(double x, double y) {
    float xf = float(WINDOW_WIDTH) / 2.0f - float(LX) / 2.0f * factor + float(x) * factor;
    float xy = float(WINDOW_HEIGHT) / 2.0f + float(LY) / 2.0f * factor - float(y) * factor;
    return {xf, xy};
}

std::ostream & operator << (std::ostream & os, sf::Vector2f const & v) {
    os << "{" << v.x << ", " << v.y << "}";
    return os;
}

#define GRID_COLOR (sf::Color(0x80, 0x80, 0x80))

void drawGrid(sf::RenderWindow & window) {
    sf::VertexArray line(sf::LinesStrip, 2);
    line[0].color = GRID_COLOR;
    line[1].color = GRID_COLOR;
    for (int i = 0; i <= NX; i ++) {
        line[0].position = phys2sf(double(i)*DX, 0);
        line[1].position = phys2sf(double(i)*DX, LY);
        window.draw(line);
    }
    for (int j = 0; j <= NY; j ++) {
        line[0].position = phys2sf(0, double(j)*DY);
        line[1].position = phys2sf(LX, double(j)*DY);
        window.draw(line);
    }
}

class mat2 {
public:
    mat2(mat2 const & m) : x(m.x), y(m.y) {}
    mat2(double x, double y) : x(x), y(y) {}
    mat2() : x(0.0), y(0.0) {}

    mat2 const & operator += (mat2 const & m) {
        x += m.x;
        y += m.y;
        return *this;
    }
    mat2 operator * (double v) const {
        return {x*v, y*v};
    }
    mat2 operator / (double v) const {
        return {x/v, y/v};
    }
    mat2 operator + (mat2 const & m) const {
        return {x+m.x, y+m.y};
    }
    mat2 operator - (mat2 const & m) const {
        return {x-m.x, y-m.y};
    }
    mat2 operator - () const {
        return {-x, -y};
    }
    [[nodiscard]] double abs() const {
        return sqrt(x*x + y*y);
    }
    [[nodiscard]] double getX() const {
        return x;
    }
    [[nodiscard]] double getY() const {
        return y;
    }

private:
    double x, y;
    friend double dot(mat2 const & m1, mat2 const & m2);
    friend std::ostream & operator << (std::ostream & os, mat2 const & m);
};

std::ostream & operator << (std::ostream & os, mat2 const & m) {
    os << "{" << m.x << ", " << m.y << "}";
    return os;
}

double dot(mat2 const & m1, mat2 const & m2) {
    return m1.x*m2.x + m1.y*m2.y;
}

class Particle : public sf::CircleShape {
public:
    Particle(mat2 const & X, mat2 const & V, Particle ** world, sf::Color const & color) : state{X, V}, buffer{X, V}, world(world) {
        setRadius(sf_diam/2.0f);
        setOrigin(sf_diam/2.0f, sf_diam/2.0f);
        setFillColor(color);
        setPosition(phys2sf(state.X.getX(), state.X.getY()));
    }

    void step(double dt) {
        if (!committed) {
            std::cerr << "Stepping with uncommitted previous step" << std::endl;
            exit(EXIT_FAILURE);
        }

        buffer.X += state.V*dt;
        buffer.V += (f_c() + f_b())/M*dt;

        if ((state.X.getY() <= D/2.0 && buffer.V.getY() < 0) || (state.X.getY() >= LY-D/2.0 && buffer.V.getY() > 0))
            buffer.V = {buffer.V.getX(), -buffer.V.getY()*0.1};
        if ((state.X.getX() <= D/2.0 && buffer.V.getX() < 0) || (state.X.getX() >= LX-D/2.0 && buffer.V.getX() > 0))
            buffer.V = {-buffer.V.getX()*0.1, buffer.V.getY()};

        committed = false;
    }

    void commit() {
        if (committed) {
            std::cerr << "No changes to commit" << std::endl;
            exit(EXIT_FAILURE);
        }

        state.X = buffer.X;
        state.V = buffer.V;

        committed = true;
        setPosition(phys2sf(state.X.getX(), state.X.getY()));
    }

private:
    double d(Particle const * p) const {
        return (state.X - p->state.X).abs();
    }

    mat2 n(Particle const * p) const {
        return (p->state.X - state.X)/d(p);
    }

    double delta(Particle const * p) const {
        return D - d(p);
    }

    mat2 v_n(Particle const * p) const {
        mat2 norm = n(p);
        return norm * dot(state.V - p->state.V, norm);
    }

    mat2 f_c();

    static mat2 f_b() {
        return {0.0, -g};
    }

    struct {
        mat2 X, V;
    } state;
    struct {
        mat2 X, V;
    } buffer;
    bool committed = true;
    static float sf_diam;

    Particle ** world;
    Particle * next = nullptr;
    int chunk = -1;

    friend void get_chunk(Particle * p, int & i, int & j);
    friend void insert_particle(Particle * p, Particle * world[NX][NY]);
    friend void update_particle(Particle * p, Particle * world[NX][NY]);
    friend void print_world(std::ostream & os, Particle * world[NX][NY]);

    friend std::ostream & operator << (std::ostream & os, Particle const & p);
};

std::ostream & operator << (std::ostream & os, Particle const & p) {
    os << "### Particle " << &p << " ###" << std::endl;
    os << "state: " << p.state.X << " " << p.state.V << std::endl;
    os << "buffer: " << p.buffer.X << " " << p.state.V << std::endl;
    os << "chunk: " << p.chunk << "/" << NY*NX-1 << std::endl;
    os << "committed: " << (p.committed ? "true" : "false") << std::endl;
    os << "next: " << p.next;
    return os;
}

float Particle::sf_diam = float(D) * factor;

void get_chunk(Particle * p, int & i, int & j) {
    i = int(p->state.X.getX()/DX);
    j = int(p->state.X.getY()/DY);
}

void insert_particle(Particle * p, Particle * world[NX][NY]) {
    int i, j;
    get_chunk(p, i, j);
    p->chunk = i*NY + j;
    p->next = world[i][j];
    world[i][j] = p;
}

void print_world(std::ostream & os, Particle * world[NX][NY]) {
    Particle ** ptr = &world[0][0];
    for (int i = 0; i < NX*NY; i ++) {
        os << "Chunk " << i << ": ";
        if (ptr[i] == nullptr) {
            os << "---" << std::endl;
        } else {
            Particle * node = ptr[i];
            while (node != nullptr) {
                os << node->state.X << " ";
                node = node->next;
            }
            os << std::endl;
        }
    }
}

void update_particle(Particle * p, Particle * world[NX][NY]) {
    int i, j;
    get_chunk(p, i, j);
    int new_chunk = i*NY + j;

    if (p->chunk < 0 || p->chunk >= NX*NY) {
        std::cerr << "Error updating chunks: untracked instance" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (new_chunk == p->chunk)
        return;

    Particle ** ptr = &world[0][0];
    Particle * p1 = ptr[p->chunk];
    if (p1 == p) {
        ptr[p->chunk] = p1->next;
    } else if (p1->next == nullptr && p1 != p) {
        std::cerr << "Error updating chunks: nullptr reached (1)" << std::endl;
        exit(EXIT_FAILURE);
    } else {
        Particle * p2 = p1->next;
        while (p2 != nullptr && p2 != p) {
            p1 = p1->next;
            p2 = p2->next;
        }
        if (p2 == nullptr) {
            std::cerr << *p << std::endl;
            print_world(std::cerr, world);
            std::cerr << "Error updating chunks: nullptr reached (2)" << std::endl;
            exit(EXIT_FAILURE);
        }
        p1->next = p2->next;
    }

    p->next = world[i][j];
    p->chunk = new_chunk;
    world[i][j] = p;
}

#define COND_CHUNK(__COND__, __I__, __J__) if (__COND__) { *buff = world[__I__][__J__]; buff ++; }

void get_neighbours(Particle * p, Particle ** buff, Particle * world[NX][NY]) {
    int i, j;
    get_chunk(p, i, j);

    *buff = world[i][j];
    buff ++;
    COND_CHUNK(i > 0, i-1, j);
    COND_CHUNK(i < NX-1, i+1, j);
    COND_CHUNK(j > 0, i, j-1);
    COND_CHUNK(j < NY-1, i, j+1);
    COND_CHUNK(i > 0 && j > 0, i-1, j-1);
    COND_CHUNK(i > 0 && j < NY-1, i-1, j+1);
    COND_CHUNK(i < NX-1 && j > 0, i+1, j-1);
    COND_CHUNK(i < NX-1 && j < NY-1, i+1, j+1);
}

#define RESTORE_ARR(__PTR__) (*reinterpret_cast<Particle *(*)[NX][NY]>(__PTR__))

mat2 Particle::f_c() {
    mat2 res;
    Particle * chunks[16] = {nullptr};
    get_neighbours(this, chunks, RESTORE_ARR(world));

    for (Particle * c : chunks) {
        if (c == nullptr)
            continue;
        while (c != nullptr) {
            if (c != this) {
                double del = delta(c);
                if (del > 0.0) {
                    res += -n(c)*del*kc - v_n(c)*eta;
                }
            }
            c = c->next;
        }
    }

    return res;
}

int main() {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(D/2.0, LX-D/2.0);
    std::uniform_int_distribution<int> dist_int(0, 0xFF);

    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Window", sf::Style::Default, sf::ContextSettings(0, 0, 4));

    std::vector<Particle *> particles;
    Particle * world[NX][NY] = {nullptr};
    auto p = new Particle({dist(mt), H0-D/2.0}, {0.0, 0.0}, &world[0][0], sf::Color(dist_int(mt), dist_int(mt), dist_int(mt)));
    particles.emplace_back(p);
    insert_particle(p, world);

    unsigned count = 0;

    while (window.isOpen()) {
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        count ++;
        if (count == 1'000) {
            count = 0;
            auto pp = new Particle({dist(mt), H0-D/2.0}, {0.0, 0.0}, &world[0][0], sf::Color(dist_int(mt), dist_int(mt), dist_int(mt)));
            particles.emplace_back(pp);
            insert_particle(pp, world);
        }

        window.clear(sf::Color::White);
        drawGrid(window);
        for (auto pp : particles) {
            window.draw(*pp);
            pp->step(0.00012);
        }
        for (auto pp : particles) {
            pp->commit();
            update_particle(pp, world);
        }
        window.display();
    }

    for (auto pp : particles) {
        delete pp;
    }

    return 0;
}
