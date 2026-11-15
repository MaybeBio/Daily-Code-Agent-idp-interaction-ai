public:
    std::string name;
    std::vector<double> argv;

    void init(std::string s, double h_0, double t_0);
    double operator()(double time);
    void print() const;
};
