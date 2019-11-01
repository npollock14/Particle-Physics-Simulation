import java.util.ArrayList;

public class Particle {
	public Position pos;
	public ArrayList<Position> positions = new ArrayList<Position>();
	public String type;
	public Position velocity = new Position(0,0,0);
	public ArrayList<Vector> forces = new ArrayList<Vector>();
	
	public Particle(Position Tpos) {
		pos = Tpos;
	}
	public Particle(Position Tpos, String Ttype) {
		pos = Tpos;
		type = Ttype;
	}
	public void setPos(Position pos) {
		this.pos = pos;
	}
	public void setType(String type) {
		this.type = type;
	}
	public void recordPos() {
		positions.add(pos);
		
	}
	public void setVelocity(Position velocity) {
		this.velocity = velocity;
	}
	public void move(double deltax, double deltay, double deltaz) {
		pos = new Position(pos.x + deltax, pos.y + deltay, pos.z + deltaz);
	}
	public void moveVec(Position vec) {
		pos = new Position(pos.x + vec.x, pos.y + vec.y, pos.z + vec.z);
	}
	
}