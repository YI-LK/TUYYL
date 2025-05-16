import matplotlib.pyplot as plt

def read_and_plot(filename="output.txt"):
    fig, ax = plt.subplots()
    points = []
    circles = []
    lines = []
    segments = []  # 新增线段列表

    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if parts[0] == 'point':
                x, y = float(parts[1]), float(parts[2])
                points.append((x, y))
            elif parts[0] == 'circle':
                x, y, r = float(parts[1]), float(parts[2]), float(parts[3])
                circles.append((x, y, r))
            elif parts[0] == 'line':
                x, y, dx, dy = map(float, parts[1:])
                length = 10
                norm = (dx**2 + dy**2)**0.5
                if norm > 0:
                    dx_norm = dx / norm * length
                    dy_norm = dy / norm * length
                else:
                    dx_norm, dy_norm = 0, 0
                x1, y1 = x - dx_norm, y - dy_norm
                x2, y2 = x + dx_norm, y + dy_norm
                lines.append(((x1, y1), (x2, y2)))
            elif parts[0] == 'segment':  # 新增线段类型
                x1, y1, x2, y2 = map(float, parts[1:])
                segments.append(((x1, y1), (x2, y2)))

    # 绘制点
    for x, y in points:
        ax.plot(x, y, 'bo')

    # 绘制圆
    for x, y, r in circles:
        circle = plt.Circle((x, y), r, color='red', fill=False)
        ax.add_patch(circle)

    # 绘制直线
    for (p1, p2) in lines:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'g-', linewidth=1)

    # 绘制线段
    for (p1, p2) in segments:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'orange', linewidth=2)  # 橙色线段

    # 设置坐标轴等比例
    ax.set_aspect('equal')
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Points, Circles, Lines and Segments")
    plt.grid(True)
    plt.savefig("output.png")
    plt.show()

if __name__ == "__main__":
    read_and_plot()