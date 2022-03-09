#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkImageHistogram.h>
#include <vtkMetaImageReader.h>
#include <vtkMarchingCubes.h>
#include <vtkConeSource.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkChartHistogram2D.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageReader2.h>
#include <vtkPlane.h>
#include <vtkBarChartActor.h> //Actor especializado para gráficos de barra
#include <vtkProperty2D.h> //vtkProperty2D contiene propiedades usadas para render imágenes y anotaciones en dos dimensiones
#include <vtkFieldData.h> //vtkFieldData representa y manipula campos de datos
#include <vtkLegendBoxActor.h> //Actor especializado para rótulos
#include <vtkStripper.h> //Es un filtro que genera tiras triangulares y/o poly-lines desde una entrada de polígonos
#include <vtkClipPolyData.h> //Es un filtro que recorta datos poligonales utilizando cualquier subclase de vtkImplicitFunction
#include <vtkImageData.h> //Leer información de imágenes
#include <vtkIntArray.h> //Arreglos
#include <vtkImageExtractComponents.h> //Extraer componentes de colores de una imágen
#include <vtkImageAccumulate.h> //Este filtro divide el espacio de componentes en contenedores discretos
#include <vtkImplicitPlaneWidget2.h> //Definir un palno interactivo infinito 
#include <vtkImplicitPlaneRepresentation.h> //Es una representación de vtkImplicitPlaneWidget2
#include <vtkFlyingEdges3D.h> //Generar isosuperficie a partir de datos de una imágen 3D
#include <vtkMatrix4x4.h> //Clase que representa una matriz de 4x4
#include <vtkImageReslice.h> //Filtro de geometría de imágen
#include <vtkImageMapToColors.h> //Toma una imagen de entrada de cualquier tipo escalar válido
#include <vtkLookupTable.h>// Objeto que es usado por los objetos mapper hacia valores escalares dentro de RGBA
#include <vtkImageActor.h>//Actor especializado para imágenes
#include <vtkInteractorStyleImage.h>//Permite usar interactivamente de la cámara (rotate, pam, zoom, etc.)

#define vtkSPtr vtkSmartPointer
#define vtkSPtrNew(Var, Type) vtkSPtr<Type> Var = vtkSPtr<Type>::New();
using namespace std;

int main() {
	vtkSPtrNew(colors, vtkNamedColors);
	vtkSPtrNew(reader, vtkMetaImageReader);	

	reader->SetFileName("../DataSources/FullHead.mhd");
	reader->Update();
// top
	// Leer archivo y mostrar
	vtkSPtrNew(boneExtractor, vtkMarchingCubes);
	boneExtractor->SetInputConnection(reader->GetOutputPort());
	boneExtractor->SetNumberOfContours(1);
	boneExtractor->SetValue(1, 1150);

	vtkSPtrNew(boneStripper, vtkStripper);
	boneStripper->SetInputConnection(boneExtractor->GetOutputPort());

	vtkSPtrNew(mapper, vtkPolyDataMapper);
	mapper->SetInputConnection(boneStripper->GetOutputPort());
	mapper->ScalarVisibilityOff(); //Quita el color escalar

	vtkSPtrNew(actor2, vtkActor);
	actor2->SetMapper(mapper);	
	actor2->GetProperty()->SetDiffuseColor(colors->GetColor3d("Ivory").GetData());	

	vtkSPtrNew(renderer2, vtkRenderer);	
	renderer2->AddActor(actor2);
	renderer2->SetBackground(colors->GetColor3d("LightGrey").GetData());
	renderer2->GetActiveCamera()->SetPosition(12, -6, 0);
	renderer2->SetViewport(0, 0.5, 1, 1); //(xmin,ymin,xmax,ymax)

// bottom right
	// Leer archivo
	// Mostrar histograma
	vtkNew<vtkImageReader2Factory> readerFactory;
	vtkSmartPointer<vtkImageReader2> reader1;
	reader1.TakeReference(readerFactory->CreateImageReader2("../DataSources/FullHead.mhd"));
	reader1->SetFileName("../DataSources/FullHead.mhd");
	reader1->Update();

	int numComponents = reader1->GetOutput()->GetNumberOfScalarComponents();
	std::cout << "Number of components: " << numComponents << std::endl;
	
	double colores[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };

	vtkNew<vtkIntArray> redFrequencies;
	//vtkNew<vtkIntArray> greenFrequencies;
	//vtkNew<vtkIntArray> blueFrequencies;

	//Ideal extraer los tres componentes, sin embargo, para el taller la imágen solo tiene un componente
	vtkNew<vtkImageExtractComponents> extract;
	extract->SetInputConnection(reader1->GetOutputPort());
	extract->SetComponents(0);
	extract->Update();

	vtkSPtrNew(histogram, vtkImageAccumulate);
	histogram->SetInputConnection(extract->GetOutputPort());
	histogram->SetComponentExtent(1, 55, 0, 0, 0, 0);
	histogram->SetComponentOrigin(0, 0, 0);
	histogram->SetComponentSpacing(1, 0, 0);
	histogram->SetIgnoreZero(0);
	histogram->Update();

	vtkIntArray* currentArray = 0;
	currentArray = redFrequencies;

	currentArray->SetNumberOfComponents(1);
	currentArray->SetNumberOfTuples(54);

	vtkIdType* output =
		static_cast<vtkIdType*>(histogram->GetOutput()->GetScalarPointer());

	for (int j = 0; j < 54; ++j)
	{
		currentArray->SetTuple1(j, *output++);
	}

	vtkNew<vtkDataObject> dataObject;
	dataObject->GetFieldData()->AddArray(redFrequencies);

	vtkNew<vtkBarChartActor> barChart;
	barChart->SetInput(dataObject);
	barChart->SetTitle("Histograma");
	barChart->GetPositionCoordinate()->SetValue(0.05, 0.05, 0.0);
	barChart->GetPosition2Coordinate()->SetValue(0.95, 0.85, 0.0);
	barChart->GetProperty()->SetColor(1, 1, 1);

	barChart->GetLegendActor()->SetNumberOfEntries(
		dataObject->GetFieldData()->GetArray(0)->GetNumberOfTuples());
	barChart->LegendVisibilityOff();
	barChart->LabelVisibilityOff();

	int count = 0;
	for (int i = 0; i < 54; ++i)
	{
		for (int j = 0; j < numComponents; ++j)
		{
			barChart->SetBarColor(count++, colores[j]);
		}
	}

	//Visualizar el histograma(s)
	vtkSPtrNew(renderer1, vtkRenderer);
	renderer1->AddActor (barChart);
	renderer1->SetBackground (colors->GetColor3d("LightSteelBlue").GetData());
	renderer1->GetActiveCamera()->SetPosition (0,0,6);
	renderer1->GetActiveCamera()->SetParallelProjection (true);
	renderer1->SetViewport (0.5, 0.0, 1, 0.5); //(xmin,ymin,xmax,ymax)

// bottom left// 
	// Leer archivo y mostrar sección (capa)
	int extent[6];
	double spacing[3];
	double origin[3];

	reader->GetOutput()->GetExtent(extent);
	reader->GetOutput()->GetSpacing(spacing);
	reader->GetOutput()->GetOrigin(origin);

	double center[3];
	center[0] = origin[0] + spacing[0] * 0.5 * (extent[0] + extent[1]);
	center[1] = origin[1] + spacing[1] * 0.5 * (extent[2] + extent[3]);
	center[2] = origin[2] + spacing[2] * 0.5 * (extent[4] + extent[5]);

	static double axialElements[16] = {
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};

	vtkSmartPointer<vtkMatrix4x4> resliceAxes =
		vtkSmartPointer<vtkMatrix4x4>::New();
	resliceAxes->DeepCopy(axialElements);
	
	resliceAxes->SetElement(0, 3, center[0]);
	resliceAxes->SetElement(1, 3, center[1]);
	resliceAxes->SetElement(2, 3, center[2]);

	vtkSmartPointer<vtkImageReslice> reslice =
		vtkSmartPointer<vtkImageReslice>::New();
	reslice->SetInputConnection(reader->GetOutputPort());
	reslice->SetOutputDimensionality(2);
	reslice->SetResliceAxes(resliceAxes);
	reslice->SetInterpolationModeToLinear();


	vtkSmartPointer<vtkLookupTable> colorTable =
		vtkSmartPointer<vtkLookupTable>::New();
	colorTable->SetRange(0, 1000);
	colorTable->SetValueRange(0.0, 1.0);
	colorTable->SetSaturationRange(0.0, 0.0);
	colorTable->SetRampToLinear();
	colorTable->Build();

	vtkSmartPointer<vtkImageMapToColors> colorMap =
		vtkSmartPointer<vtkImageMapToColors>::New();
	colorMap->SetLookupTable(colorTable);
	colorMap->SetInputConnection(reslice->GetOutputPort());
	colorMap->Update();colorMap->Update();

	vtkSmartPointer<vtkImageActor> imgActor =
		vtkSmartPointer<vtkImageActor>::New();
	imgActor->SetInputData(colorMap->GetOutput());	

	vtkSPtrNew(renderer0, vtkRenderer);	
	renderer0->AddActor(imgActor);
	renderer0->SetBackground(colors->GetColor3d("LightYellow").GetData());
	renderer0->GetActiveCamera()->SetPosition(0, 0, 6);
	renderer0->GetActiveCamera()->SetParallelProjection(true);
	renderer0->SetViewport(0.0, 0.0, 0.5, 0.5); //(xmin,ymin,xmax,ymax)

// Ventana completa
	renderer0->ResetCamera();
	renderer1->ResetCamera();

	renderer2->GetActiveCamera()->Azimuth(30);
	renderer2->GetActiveCamera()->Elevation(30);
	renderer2->ResetCamera();
	renderer2->GetActiveCamera()->Zoom(0.75);	

	vtkSPtrNew(renderWindow, vtkRenderWindow);
	renderWindow->Render();

	renderWindow->AddRenderer(renderer2);
	renderWindow->AddRenderer (renderer0);
	renderWindow->AddRenderer (renderer1);	
	renderWindow->SetSize (1600,800);

	vtkSPtrNew(renderWindowInteractor, vtkRenderWindowInteractor);

	vtkNew<vtkImplicitPlaneRepresentation> rep;
	rep->SetPlaceFactor(1.25); // Tamaño de la caja widget
	rep->PlaceWidget(actor2->GetBounds());
	vtkNew<vtkPlane> plane;
	rep->SetNormal(plane->GetNormal());

	vtkNew<vtkImplicitPlaneWidget2> planeWidget;
	planeWidget->SetInteractor(renderWindowInteractor);
	planeWidget->SetRepresentation(rep);
	planeWidget->SetCurrentRenderer(renderer2);

	vtkSmartPointer<vtkInteractorStyleImage> imagestyle =
		vtkSmartPointer<vtkInteractorStyleImage>::New();
	renderWindowInteractor->SetInteractorStyle(imagestyle);
	renderWindowInteractor->SetRenderWindow (renderWindow);

	renderWindow->Render ();

	planeWidget->On();
	renderWindowInteractor->Start();
	return 0;
}