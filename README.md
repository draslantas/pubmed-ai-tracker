# Kişisel PubMed AI Asistanı 🧬

Bu program, PubMed veritabanını tarayarak "Artificial Intelligence" konulu yeni makaleleri bulur, kategorize eder ve size sunar.

## 📦 Kurulum
Bu klasörde bir terminal açın ve gerekli kütüphaneleri yükleyin:
```bash
pip3 install streamlit schedule biopython
```

## 🚀 Nasıl Çalıştırılır?

### 1. Arayüzü (Dashboard) Açmak İçin
Mavi/Yeşil kartların olduğu ekranı açmak için şu komutu kullanın:
```bash
streamlit run app.py
```

### 2. Otomatik Zamanlayıcıyı (Sabah 08:00) Başlatmak İçin
Arka planda çalışıp her sabah tarama yapması için şu komutu kullanın:
```bash
python3 scheduler.py
```

## ✨ Özellikler
- **Otomatik Takip:** Her gün son 24 saatteki makaleleri bulur.
- **Akıllı Etiketleme:** Makaleleri branş (Beyin Cerrahisi, Kardiyoloji vb.) ve AI Yöntemi (Deep Learning, LLM vb.) olarak etiketler.
- **Okuma Listesi:** Okumadıklarınız **Mavi**, okuduklarınız **Yeşil** görünür.
