import pandas as pd

def eliminar_columnas_por_indice(input_csv, output_csv, indices_a_eliminar):
    df = pd.read_csv(input_csv)
    columnas_a_eliminar = [df.columns[i] for i in indices_a_eliminar]
    df = df.drop(columns=columnas_a_eliminar, errors='ignore')
    df.to_csv(output_csv, index=False)

    print(f"Archivo guardado en: {output_csv}")

if __name__ == "__main__":
    input_csv = "abalone.csv"
    output_csv = "abalone_real.csv"
    indices_a_eliminar = [0,7,8]  
    eliminar_columnas_por_indice(input_csv, output_csv, indices_a_eliminar)
