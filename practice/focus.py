from photometry import photometry
from simulation import simulation

class Image:
    def __init__(self, focus):
        self.focus = focus

        true_fwhm = 0.016*(focus-360)**2+6.0
        simulation(fwhm=true_fwhm)

        file_name = str(true_fwhm)
        file_name = file_name.replace('.', '_')
        # self.fwhm = photometry(f'images/focus_{file_name}.fits')
        


def focus_finder(image1, image2, seen):

    print(f"image1: focus:{image1.focus}, FWHM:{image1.fwhm}")
    print(f"image2: focus:{image2.focus}, FWHM:{image2.fwhm} \n")

    if image1.fwhm > image2.fwhm:
        focus3 = image2.focus + (image2.focus - image1.focus)
        if focus3 in seen:
            return image2.focus
        seen.add(focus3)
        image3 = Image(focus3)
        return focus_finder(image2, image3, seen)
    elif image1.fwhm < image2.fwhm:
        focus3 = image1.focus - (image2.focus - image1.focus)
        if focus3 in seen:
            return image1.focus
        seen.add(focus3)
        image3 = Image(focus3)
        return focus_finder(image3, image1, seen)
    else:
        return image1.focus

def curve_finder(image1, image2, seen, direction=None):
    
    if image1.fwhm > image2.fwhm:
        if direction is None:
            direction = 'right'
        if direction == 'left':
            return curve_helper(image1, image2, seen)
        focus3 = image2.focus + (image2.focus - image1.focus)
        image3 = Image(focus3)
        seen.add(image3)
        return curve_finder(image2, image3, seen, direction)
    elif image1.fwhm < image2.fwhm:
        if direction is None:
            direction = 'left'
        if direction == 'right':
            return curve_helper(image1, image2, seen)
        focus3 = image1.focus - (image2.focus - image1.focus)
        image3 = Image(focus3)
        seen.add(image3)
        return curve_finder(image3, image1, seen, direction)
    
def curve_helper(image1, image2, seen, iterations=2):
    if iterations > 0:
        if image1.fwhm > image2.fwhm:
            focus3 = image1.focus - (image2.focus - image1.focus)
            image3 = Image(focus3)
            seen.add(image3)
            return curve_helper(image3, image1, seen, iterations-1)
        elif image1.fwhm < image2.fwhm:
            focus3 = image2.focus + (image2.focus - image1.focus)
            image3 = Image(focus3)
            seen.add(image3)
            return curve_helper(image2, image3, seen, iterations-1)
    else:
        return seen

    
import argparse
import quadratic

def main():
    parser = argparse.ArgumentParser(description="Find optimal focus for images based on FWHM.")
    parser.add_argument('--image1', type=int, default=380, help='Focus value for the first image')
    parser.add_argument('--image2', type=int, default=385, help='Focus value for the second image')
    args = parser.parse_args()
    # image1 = Image(335)
    # image2 = Image(340)
    # image1 = Image(380)
    # image2 = Image(385)

    image1 = Image(args.image1)
    image2 = Image(args.image2)

    seen = set()
    seen.add(image1.focus)
    seen.add(image2.focus)

    focus = focus_finder(image1, image2, seen)
    print(f"Optimal focus: {focus}\n")

    seen2 = set()
    seen2.add(image1)
    seen2.add(image2)
    curve = curve_finder(image1, image2, seen2)
    curve = sorted(curve, key=lambda img: img.focus)
    x_values = []
    y_values = []
    print("Curve found with the following focus values:")
    for img in curve:
        print(f"Focus: {img.focus}, FWHM: {img.fwhm}")
        x_values.append(img.focus)
        y_values.append(img.fwhm)

    a, b, c = quadratic.fit_quadratic(x_values, y_values)
    x_vertex, y_vertex = quadratic.vertex(a, b, c)
    print(f"Optimal focus: {x_vertex}, FWHM: {y_vertex}")


if __name__ == "__main__":
    main()
