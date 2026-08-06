import numpy as np
import random


def momentum_angles(px, py, pz):
    """
    Calculate polar and azimuthal angles from the momentum vector.
    """
    p = np.sqrt(px**2 + py**2 + pz**2)

    theta = np.arccos(pz / p)
    phi = np.arctan2(py, px)

    return theta, phi


def find_max_weight(filename):
    """
    First pass over the LHE file to determine the maximum event weight.
    """

    w_max = 0.0

    with open(filename, "r") as f:

        for line in f:

            if "<event" in line:

                header = next(f).split()

                weight = float(header[2])

                if abs(weight) > w_max:
                    w_max = abs(weight)

    return w_max



def parse_lhe(filename, outputfile, w_max):
    """
    Parse LHE events and perform accept/reject unweighting.
    """

    accepted = 0
    total = 0

    with open(filename, "r") as f, open(outputfile, "w") as out:

        for line in f:

            if "<event" in line:

                total += 1

                # Event header:
                # NUP IDPRUP XWGTUP ...
                header = next(f).split()

                n_particles = int(header[0])
                weight = float(header[2])


                # Accept/reject step
                probability = weight / w_max

                if random.random() > probability:

                    # Skip this event
                    for _ in range(n_particles):
                        next(f)

                    continue


                muon = None
                proton = None
                photon = None


                for _ in range(n_particles):

                    data = next(f).split()

                    pid = int(data[0])

                    # Four momentum
                    px = float(data[6])
                    py = float(data[7])
                    pz = float(data[8])

                    E = float(data[9])
                    mass = float(data[10])


                    theta, phi = momentum_angles(px, py, pz)


                    # Outgoing muon
                    if pid == 13:
                        muon = [E, theta, phi]


                    # Recoil proton
                    elif pid == 2212:
                        proton = [E, theta, phi]


                    # Bremsstrahlung photon
                    elif pid == 0 and abs(mass) < 1e-10:
                        photon = [E, theta, phi]


                if muon is not None and proton is not None:

                    if photon is None:
                        photon = [0.0, 0.0, 0.0]


                    values = (
                        muon
                        + proton
                        + photon
                    )


                    out.write(
                        " ".join(f"{x:.6f}" for x in values)
                        + "\n"
                    )

                    accepted += 1


    print(f"Total events: {total}")
    print(f"Accepted events: {accepted}")
    print(f"Acceptance rate: {accepted/total:.3f}")



if __name__ == "__main__":

    inputfile = "events.lhe"
    outputfile = "events.txt"


    print("Finding maximum weight...")
    w_max = find_max_weight(inputfile)

    print(f"Maximum weight: {w_max}")


    print("Parsing events...")
    parse_lhe(
        inputfile,
        outputfile,
        w_max
    )