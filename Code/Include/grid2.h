#include <vector>

template<typename T>
class Grid2
{
public:
	inline Grid2() { }
	inline Grid2(int nx, int ny) : nx(nx), ny(ny)
	{
		data.resize(nx * ny);
	}
	inline Grid2(int nx, int ny, const T& defaultVal) : nx(nx), ny(ny)
	{
		data.resize(nx * ny, defaultVal);
	}

	inline T& operator()(int i, int j)
	{
		return data[ToIndex1D(i, j)];
	}

	inline T& operator()(int i)
	{
		return data[i];
	}
	inline T operator()(int i, int j) const
	{
		return data[ToIndex1D(i, j)];
	}

	inline T& operator[](int i) {
		return data[i];
	}
	inline T operator[](int i) const {
		return data[i];
	}

	inline int ToIndex1D(int i, int j) const
	{
		return j * nx + i;
	}
	inline void ToIndex2D(int k, int& i, int& j) const
	{
		j = k / nx;
		i = k % nx;
	}
	inline bool Inside(int i, int j) const
	{
		return (i >= 0 && i < nx && j >= 0 && j < ny);
	}

	inline void* GetData() { return data.data(); }
	inline int Size() const { return nx * ny; }
	inline size_t ByteSize() const { return Size() * sizeof(T); }
	inline int Width() const { return nx; }
	inline int Height() const { return ny; }

protected:
	std::vector<T> data;
	int nx, ny;
};
